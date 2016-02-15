#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#
# Ref:
# J. Chem. Phys. 117, 7433
#

import time
from functools import reduce
import numpy
import pyscf.lib
from pyscf.lib import logger
from pyscf.dft import numint
from pyscf.scf import cphf
from pyscf.tddft import rks
from pyscf.tddft import rhf_grad
import pyscf.dft.vxc


#
# Given Y = 0, TDDFT gradients (XAX+XBY+YBX+YAY)^1 turn to TDA gradients (XAX)^1
#
def kernel(td_grad, (x, y), singlet=True, atmlst=None,
           max_memory=2000, verbose=logger.INFO):
    if isinstance(verbose, logger.Logger):
        log = verbose
    else:
        log = logger.Logger(td_grad.stdout, verbose)
    time0 = time.clock(), time.time()

    mol = td_grad.mol
    mf = td_grad._td._scf
    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    mo_occ = mf.mo_occ
    nao, nmo = mo_coeff.shape
    nocc = (mo_occ>0).sum()
    nvir = nmo - nocc
    xpy = (x+y).reshape(nvir,nocc)
    xmy = (x-y).reshape(nvir,nocc)
    orbv = mo_coeff[:,nocc:]
    orbo = mo_coeff[:,:nocc]

    dvv = numpy.einsum('ai,bi->ab', xpy, xpy) + numpy.einsum('ai,bi->ab', xmy, xmy)
    doo =-numpy.einsum('ai,aj->ij', xpy, xpy) - numpy.einsum('ai,aj->ij', xmy, xmy)
    dmzvop = reduce(numpy.dot, (orbv, xpy, orbo.T))
    dmzvom = reduce(numpy.dot, (orbv, xmy, orbo.T))
    dmzoo = reduce(numpy.dot, (orbo, doo, orbo.T))
    dmzoo+= reduce(numpy.dot, (orbv, dvv, orbv.T))

    mem_now = pyscf.lib.current_memory()[0]
    max_memory = max(2000, td_grad.max_memory*.9-mem_now)

    x_code, c_code = pyscf.dft.vxc.parse_xc_name(mf.xc)
    hyb = mf._numint.hybrid_coeff(x_code, spin=(mol.spin>0)+1)
    f1vo, f1oo, vxc1, k1ao = \
            _contract_xc_kernel(td_grad, x_code, c_code, dmzvop,
                                dmzoo, True, True, singlet, max_memory)

    if abs(hyb) > 1e-10:
        vj, vk = mf.get_jk(mol, (dmzoo, dmzvop+dmzvop.T, dmzvom-dmzvom.T), hermi=0)
        veff0doo = vj[0] * 2 - hyb * vk[0] + f1oo[0] + k1ao[0] * 2
        wvo = reduce(numpy.dot, (orbv.T, veff0doo, orbo)) * 2
        if singlet:
            veff = vj[1] * 2 - hyb * vk[1] + f1vo[0] * 2
        else:
            veff = -hyb * vk[1]
        veff0mop = reduce(numpy.dot, (mo_coeff.T, veff, mo_coeff))
        wvo -= numpy.einsum('ki,ai->ak', veff0mop[:nocc,:nocc], xpy) * 2
        wvo += numpy.einsum('ac,ai->ci', veff0mop[nocc:,nocc:], xpy) * 2
        veff = -hyb * vk[2]
        veff0mom = reduce(numpy.dot, (mo_coeff.T, veff, mo_coeff))
        wvo -= numpy.einsum('ki,ai->ak', veff0mom[:nocc,:nocc], xmy) * 2
        wvo += numpy.einsum('ac,ai->ci', veff0mom[nocc:,nocc:], xmy) * 2
    else:
        vj = mf.get_j(mol, (dmzoo, dmzvop+dmzvop.T), hermi=1)
        veff0doo = vj[0] * 2 + f1oo[0] + k1ao[0] * 2
        wvo = reduce(numpy.dot, (orbv.T, veff0doo, orbo)) * 2
        if singlet:
            veff = vj[1] * 2 + f1vo[0] * 2
        else:
            veff = f1vo[0] * 2
        veff0mop = reduce(numpy.dot, (mo_coeff.T, veff, mo_coeff))
        wvo -= numpy.einsum('ki,ai->ak', veff0mop[:nocc,:nocc], xpy) * 2
        wvo += numpy.einsum('ac,ai->ci', veff0mop[nocc:,nocc:], xpy) * 2
        veff0mom = numpy.zeros((nmo,nmo))
    def fvind(x):
# Cannot make call to ._td.get_vind because first order orbitals are solved
# through closed shell ground state CPHF.
        dm = reduce(numpy.dot, (orbv, x.reshape(nvir,nocc), orbo.T))
# Call singlet XC kernel contraction, for closed shell ground state
        vindxc = rks._contract_xc_kernel(td_grad._td, x_code, c_code,
                                         [dm+dm.T], True, max_memory)
        if abs(hyb) > 1e-10:
            vj, vk = mf.get_jk(mol, (dm+dm.T))
            veff = vj * 2 - hyb * vk + vindxc
        else:
            vj = mf.get_j(mol, (dm+dm.T))
            veff = vj * 2 + vindxc
        return reduce(numpy.dot, (orbv.T, veff, orbo)).ravel()
    z1 = cphf.solve(fvind, mo_energy, mo_occ, wvo,
                    max_cycle=td_grad.max_cycle_cphf, tol=td_grad.conv_tol)[0]
    z1 = z1.reshape(nvir,nocc)
    time1 = log.timer('Z-vector using CPHF solver', *time0)

    z1ao  = reduce(numpy.dot, (orbv, z1, orbo.T))
# Note Z-vector is always associated to singlet integrals.
    fxcz1 = _contract_xc_kernel(td_grad, x_code, c_code, z1ao, None,
                                False, False, True, max_memory)[0]
    if abs(hyb) > 1e-10:
        vj, vk = mf.get_jk(mol, z1ao, hermi=0)
        veff = vj * 2 - hyb * vk + fxcz1[0]
    else:
        vj = mf.get_j(mol, z1ao, hermi=1)
        veff = vj * 2 + fxcz1[0]

    im0 = numpy.zeros((nmo,nmo))
    im0[:nocc,:nocc] = reduce(numpy.dot, (orbo.T, veff0doo+veff, orbo))
    im0[:nocc,:nocc]+= numpy.einsum('ak,ai->ki', veff0mop[nocc:,:nocc], xpy)
    im0[:nocc,:nocc]+= numpy.einsum('ak,ai->ki', veff0mom[nocc:,:nocc], xmy)
    im0[nocc:,nocc:] = numpy.einsum('ci,ai->ac', veff0mop[nocc:,:nocc], xpy)
    im0[nocc:,nocc:]+= numpy.einsum('ci,ai->ac', veff0mom[nocc:,:nocc], xmy)
    im0[nocc:,:nocc] = numpy.einsum('ki,ai->ak', veff0mop[:nocc,:nocc], xpy)*2
    im0[nocc:,:nocc]+= numpy.einsum('ki,ai->ak', veff0mom[:nocc,:nocc], xmy)*2

    zeta = pyscf.lib.direct_sum('i+j->ij', mo_energy, mo_energy) * .5
    zeta[nocc:,:nocc] = mo_energy[:nocc]
    zeta[:nocc,nocc:] = mo_energy[nocc:]
    dm1 = numpy.zeros((nmo,nmo))
    dm1[:nocc,:nocc] = doo
    dm1[nocc:,nocc:] = dvv
    dm1[nocc:,:nocc] = z1
    dm1[:nocc,:nocc] += numpy.eye(nocc)*2 # for ground state
    im0 = reduce(numpy.dot, (mo_coeff, im0+zeta*dm1, mo_coeff.T))

    h1 = td_grad.get_hcore(mol)
    s1 = td_grad.get_ovlp(mol)

    dmz1doo = z1ao + dmzoo
    oo0 = reduce(numpy.dot, (orbo, orbo.T))
    if abs(hyb) > 1e-10:
        vj, vk = td_grad.get_jk(mol, (oo0, dmz1doo+dmz1doo.T, dmzvop+dmzvop.T,
                                      dmzvom-dmzvom.T))
        vj = vj.reshape(-1,3,nao,nao)
        vk = vk.reshape(-1,3,nao,nao)
        if singlet:
            veff1 = vj * 2 - hyb * vk
        else:
            veff1 = numpy.vstack((vj[:2]*2-hyb*vk[:2], -hyb*vk[2:]))
    else:
        vj = td_grad.get_j(mol, (oo0, dmz1doo+dmz1doo.T, dmzvop+dmzvop.T))
        vj = vj.reshape(-1,3,nao,nao)
        veff1 = numpy.zeros((4,3,nao,nao))
        if singlet:
            veff1[:3] = vj * 2
        else:
            veff1[:2] = vj[:2] * 2
    veff1[0] += vxc1[1:]
    veff1[1] +=(f1oo[1:] + fxcz1[1:] + k1ao[1:]*2)*2 # *2 for dmz1doo+dmz1oo.T
    veff1[2] += f1vo[1:] * 2
    time1 = log.timer('2e AO integral derivatives', *time1)

    if atmlst is None:
        atmlst = range(mol.natm)
    offsetdic = td_grad.aorange_by_atom()
    de = numpy.zeros((len(atmlst),3))
    for k, ia in enumerate(atmlst):
        shl0, shl1, p0, p1 = offsetdic[ia]

        mol.set_rinv_origin_(mol.atom_coord(ia))
        h1ao = -mol.atom_charge(ia) * mol.intor('cint1e_iprinv_sph', comp=3)
        h1ao[:,p0:p1] += h1[:,p0:p1] + veff1[0,:,p0:p1]

        # Ground state gradients
        # h1ao*2 for +c.c, oo0*2 for doubly occupied orbitals
        e1  = numpy.einsum('xpq,pq->x', h1ao, oo0) * 4

        e1 += numpy.einsum('xpq,pq->x', h1ao, dmz1doo)
        e1 += numpy.einsum('xqp,pq->x', h1ao, dmz1doo)
        e1 -= numpy.einsum('xpq,pq->x', s1[:,p0:p1], im0[p0:p1])
        e1 -= numpy.einsum('xqp,pq->x', s1[:,p0:p1], im0[:,p0:p1])

        e1 += numpy.einsum('xij,ij->x', veff1[1,:,p0:p1], oo0[p0:p1])
        e1 += numpy.einsum('xij,ij->x', veff1[2,:,p0:p1], dmzvop[p0:p1,:]) * 2
        e1 += numpy.einsum('xij,ij->x', veff1[3,:,p0:p1], dmzvom[p0:p1,:]) * 2
        e1 += numpy.einsum('xji,ij->x', veff1[2,:,p0:p1], dmzvop[:,p0:p1]) * 2
        e1 -= numpy.einsum('xji,ij->x', veff1[3,:,p0:p1], dmzvom[:,p0:p1]) * 2

        de[k] = e1

    log.timer('TDDFT nuclear gradients', *time0)
    return de

# xai, oovv in AO-representation
# Note spin-trace are applied for fxc, kxc
def _contract_xc_kernel(td_grad, x_id, c_id, xai, oovv=None, with_vxc=True,
                        with_kxc=True, singlet=True, max_memory=4000):
    mol = td_grad.mol
    mf = td_grad._scf
    ni = mf._numint
    grids = mf.grids

    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    mo_occ = mf.mo_occ
    nao, nmo = mo_coeff.shape
    nocc = (mo_occ>0).sum()
    orbv = mo_coeff[:,nocc:]
    orbo = mo_coeff[:,:nocc]

    # dmvo ~ reduce(numpy.dot, (orbv, Xai, orbo.T))
    dmvo = (xai + xai.T) * .5 # because K_{ai,bj} == K_{ai,bj}

    f1vo = numpy.zeros((4,nao,nao))
    deriv = 2
    if oovv is not None:
        f1oo = numpy.zeros((4,nao,nao))
    else:
        f1oo = None
    if with_vxc:
        v1ao = numpy.zeros((4,nao,nao))
    else:
        v1ao = None
    if with_kxc:
        k1ao = numpy.zeros((4,nao,nao))
        deriv = 3
    else:
        k1ao = None

    xctype = ni._xc_type(x_id, c_id)
    ngrids = len(grids.weights)
    BLKSIZE = numint.BLKSIZE
    max_memory = max_memory - 4*nao**2*8/1e6
    blksize = min(int(max_memory/6*1e6/8/nao/BLKSIZE)*BLKSIZE, ngrids)

    if xctype == 'LDA':
        ao_deriv = 1
        for ao, non0tab, weight, coords \
                in ni.block_loop(mol, grids, nao, ao_deriv, max_memory, ni.non0tab):
            rho = ni.eval_rho2(mol, ao[0], mo_coeff, mo_occ, non0tab, 'LDA')
            vxc, fxc, kxc = ni.eval_xc(x_id, c_id, rho, 0, deriv=deriv)[1:]

            wfxc = fxc[0] * weight * 2  # *2 for alpha+beta
            if singlet:
                rho1 = ni.eval_rho(mol, ao[0], dmvo, non0tab, 'LDA')
                aow = numpy.einsum('pi,p->pi', ao[0], wfxc*rho1)
                for k in range(4):
                    f1vo[k] += numint._dot_ao_ao(mol, ao[k], aow, nao, weight.size, non0tab)
            if oovv is not None:
                rho2 = ni.eval_rho(mol, ao[0], oovv, non0tab, 'LDA')
                aow = numpy.einsum('pi,p->pi', ao[0], wfxc*rho2)
                for k in range(4):
                    f1oo[k] += numint._dot_ao_ao(mol, ao[k], aow, nao, weight.size, non0tab)
            if with_vxc:
                aow = numpy.einsum('pi,p->pi', ao[0], vxc[0]*weight)
                for k in range(4):
                    v1ao[k] += numint._dot_ao_ao(mol, ao[k], aow, nao, weight.size, non0tab)
            if singlet and with_kxc:
                aow = numpy.einsum('pi,p->pi', ao[0], kxc[0]*weight*rho1**2)
                for k in range(4):
                    k1ao[k] += numint._dot_ao_ao(mol, ao[k], aow, nao, weight.size, non0tab)
            vxc = fxc = kxc = aow = rho = rho1 = rho2 = None
        if singlet and with_kxc:
            k1ao *= 4

    elif xctype == 'GGA':
        raise NotImplementedError('GGA')
    else:
        raise NotImplementedError('meta-GGA')

    f1vo[1:] *= -1
    if f1oo is not None: f1oo[1:] *= -1
    if v1ao is not None: v1ao[1:] *= -1
    if k1ao is not None: k1ao[1:] *= -1
    return f1vo, f1oo, v1ao, k1ao


class Gradients(rhf_grad.Gradients):
    def grad_elec(self, xy, singlet, atmlst=None):
        return kernel(self, xy, singlet, atmlst, self.max_memory, self.verbose)


if __name__ == '__main__':
    from pyscf import gto
    from pyscf import scf
    from pyscf import dft
    import pyscf.tddft
    mol = gto.Mole()
    mol.verbose = 0
    mol.output = None

    mol.atom = [
        ['H' , (0. , 0. , 1.804)],
        ['F' , (0. , 0. , 0.)], ]
    mol.unit = 'B'
    mol.basis = '631g'
    mol.build()

    mf = dft.RKS(mol)
    mf.xc = 1,0
    mf.grids.prune = False
#    mf.grids.level = 6
    mf.scf()

    td = pyscf.tddft.TDA(mf)
    td.nstates = 3
    e, z = td.kernel()
    tdg = Gradients(td)
    g1 = tdg.kernel(state=2)
    print(g1)
# [[ 0  0  -9.23916667e-02]
#  [ 0  0   9.23956206e-02]]

    td = pyscf.tddft.TDDFT(mf)
    td.nstates = 3
    e, z = td.kernel()
    tdg = Gradients(td)
    g1 = tdg.kernel(state=2)
    print(g1)
# [[ 0  0  -1.31315477e-01]
#  [ 0  0   1.31319442e-01]]

    td = pyscf.tddft.TDA(mf)
    td.nstates = 3
    td.singlet = False
    e, z = td.kernel()
    tdg = Gradients(td)
    g1 = tdg.kernel(state=2)
    print(g1)
# [[ 0  0  -2.72555279e-01]
#  [ 0  0   2.72555086e-01]]

    td = pyscf.tddft.TDDFT(mf)
    td.nstates = 3
    td.singlet = False
    e, z = td.kernel()
    tdg = Gradients(td)
    g1 = tdg.kernel(state=2)
    print(g1)
# [[ 0  0  -2.72555279e-01]
#  [ 0  0   2.72555086e-01]]
