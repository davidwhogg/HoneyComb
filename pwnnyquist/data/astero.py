import numpy as np

class astero(object):

    def __init__(self):

        # Load SDSS Astero data from table 1 - KIDs, nu_max, nu_max_err,
        # dnu, dnu_err, SDSS teffs and feh
        sdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                              skiprows=30, skip_footer=1313, invalid_raise=False,
                              usecols=(0,1,2,3,4,5,6,9,10)).T
        sKID1, snu_max1, snu_max_err1, sdnu1, sdnu_err1, steff1, steff_err1, sfeh1, \
                sfeh_err1 = sdata

        # Load IRFM Astero data from table 1 - KIDs, nu_max, nu_max_err,
        # dnu, dnu_err, IRFM teffs and feh
        idata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                              skiprows=30, skip_footer=1313, invalid_raise=False,
                              usecols=(0,1,2,3,4,7,8,9,10)).T
        iKID1, inu_max1, inu_max_err1, idnu1, idnu_err1, iteff1, iteff_err1, ifeh1, \
                ifeh_err1 = sdata

        # Load Astero data from table 2 - KIDs, nu_max, nu_max_err,
        # dnu, dnu_err, teffs and feh from bruntt
        bdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                               skiprows=576, skip_footer=1228, invalid_raise=False,
                               usecols=(0,1,2,3,4,5,6,7,8)).T
        bKID1, bnu_max1, bnu_max_err1, bdnu1, bdnu_err1, bteff1, bteff_err1, bfeh1, \
                bfeh_err1 = bdata

        # load astero data from table 4 - KID, m, m_errp, m_errm,
        # r, r_errp, r_errm, rho, rho_errp, rho_errm logg, age with SDSS teffs
        sdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                              skiprows=698, skip_footer=675, invalid_raise=False,
                              usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)).T
        sKID2, sm, sm_errp, sm_errm, sr, sr_errp, sr_errm, srho, srho_errp, \
                srho_errm, slogg, slogg_errp, slogg_errm, sage, sage_errp, \
                sage_errm = sdata

        # load astero data from table 5 - KID, m, m_errpm m_errm,
        # r, r_errp, r_errm, rho, rho_errp, rho_errm logg, age with IRFM teffs
        idata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                              skiprows=1251, skip_footer=122, invalid_raise=False,
                              usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)).T
        iKID2, im, im_errp, im_errm, ir, ir_errp, ir_errm, irho, irho_errp, \
                irho_errm, ilogg, ilogg_errp, ilogg_errm, iage, iage_errp, \
                iage_errm = idata

        # load astero data from table 6 - KID, m, m_errp, m_errm,
        # r, r_errp, r_errm, rho, rho_errp, rho_errm logg, age with Bruntt teffs
        bdata = np.genfromtxt('/Users/angusr/Python/Gyro/data/ApJtable_zeros.txt',
                              skiprows=1804, invalid_raise=False,
                              usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)).T
        bKID2, bm, bm_errp, bm_errm, br, br_errp, br_errm, brho, brho_errp, \
                brho_errm, blogg, blogg_errp, blogg_errm, bage, bage_errp, \
                bage_errm = bdata

        #trim tables 4 and 5 to be the same size as table 1
        sKID, iKID = np.zeros(len(sKID2)), np.zeros(len(iKID2))
        bKID = np.zeros(len(bKID2))
        steff, steff_err = np.zeros(len(sKID2)), np.zeros(len(sKID2))
        sfeh, sfeh_err = np.zeros(len(sKID2)), np.zeros(len(sKID2))
        snu_max, snu_max_err = np.zeros(len(sKID2)), np.zeros(len(sKID2))
        sdnu, sdnu_err = np.zeros(len(sKID2)), np.zeros(len(sKID2))
        iteff, iteff_err = np.zeros(len(iKID2)), np.zeros(len(iKID2))
        ifeh, ifeh_err = np.zeros(len(iKID2)), np.zeros(len(iKID2))
        inu_max, inu_max_err = np.zeros(len(iKID2)), np.zeros(len(iKID2))
        idnu, idnu_err = np.zeros(len(iKID2)), np.zeros(len(iKID2))
        bteff, bteff_err = np.zeros(len(bKID2)), np.zeros(len(bKID2))
        bfeh, bfeh_err = np.zeros(len(bKID2)), np.zeros(len(bKID2))
        bnu_max, bnu_max_err = np.zeros(len(bKID2)), np.zeros(len(bKID2))
        bdnu, bdnu_err = np.zeros(len(bKID2)), np.zeros(len(bKID2))

        for i, kid in enumerate(sKID2):
            l = sKID1==kid
            if sum(l):
                sKID[i] = kid
                steff[i], steff_err[i] = steff1[i], steff_err1[i]
                sfeh[i], sfeh_err[i] = sfeh1[i], sfeh_err1[i]
                snu_max[i], snu_max_err[i] = snu_max1[i], snu_max_err1[i]
                sdnu[i], sdnu_err[i] = sdnu1[i], sdnu_err1[i]

        for i, kid in enumerate(iKID2):
            l = iKID1==kid
            if sum(l):
                iKID[i] = kid
                iteff[i], iteff_err[i] = iteff1[i], iteff_err1[i]
                ifeh[i], ifeh_err[i] = ifeh1[i], ifeh_err1[i]
                inu_max[i], inu_max_err[i] = inu_max1[i], inu_max_err1[i]
                idnu[i], idnu_err[i] = idnu1[i], idnu_err1[i]

        for i, kid in enumerate(bKID2):
            l = bKID1==kid
            if sum(l):
                bKID[i] = kid
                bteff[i], bteff_err[i] = bteff1[i], bteff_err1[i]
                bfeh[i], bfeh_err[i] = bfeh1[i], bfeh_err1[i]
                bnu_max[i], bnu_max_err[i] = bnu_max1[i], bnu_max_err1[i]
                bdnu[i], bdnu_err[i] = bdnu1[i], bdnu_err1[i]

        self.iKID = iKID
        self.iteff = iteff
        self.iteff_err = iteff_err
        self.ifeh = ifeh
        self.ifeh_err = ifeh_err
        self.inu_max = inu_max
        self.inu_max_err = inu_max_err
        self.idnu = idnu
        self.idnu_err = idnu_err
        self.ilogg = ilogg
        self.ilogg_errp = ilogg_errp
        self.ilogg_errm = ilogg_errm
        self.iage = iage
        self.iage_errp = iage_errp
        self.iage_errm = iage_errm
        self.im = im
        self.im_errp = im_errp
        self.im_errm = im_errm
        self.ir = ir
        self.ir_errp = ir_errp
        self.ir_errm = ir_errm
        self.irho = irho
        self.irho_errp = irho_errp
        self.irho_errm = irho_errm

        self.sKID = sKID
        self.steff = steff
        self.steff_err = steff_err
        self.sfeh = sfeh
        self.sfeh_err = sfeh_err
        self.snu_max = snu_max
        self.snu_max_err = snu_max_err
        self.sdnu = sdnu
        self.sdnu_err = sdnu_err
        self.slogg = slogg
        self.slogg_errp = slogg_errp
        self.slogg_errm = slogg_errm
        self.sage = sage
        self.sage_errp = sage_errp
        self.sage_errm = sage_errm
        self.sm = sm
        self.sm_errp = sm_errp
        self.sm_errm = sm_errm
        self.sr = sr
        self.sr_errp = sr_errp
        self.sr_errm = sr_errm
        self.srho = srho
        self.srho_errp = srho_errp
        self.srho_errm = srho_errm

        self.bKID = bKID2
        self.bteff = bteff
        self.bteff_err = bteff_err
        self.bfeh = bfeh
        self.bfeh_err = bfeh_err
        self.bnu_max = bnu_max
        self.bnu_max_err = bnu_max_err
        self.bdnu = bdnu
        self.bdnu_err = bdnu_err
        self.blogg = blogg
        self.blogg_errp = blogg_errp
        self.blogg_errm = blogg_errm
        self.bage = bage
        self.bage_errp = bage_errp
        self.bage_errm = bage_errm
        self.bm = bm
        self.bm_errp = bm_errp
        self.bm_errm = bm_errm
        self.br = br
        self.br_errp = br_errp
        self.br_errm = br_errm
        self.brho = brho
        self.brho_errp = brho_errp
        self.brho_errm = brho_errm

data = astero()
