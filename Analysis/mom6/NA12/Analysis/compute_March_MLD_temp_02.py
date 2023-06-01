import gsw

zl = np.copy(df.zl)
pr = gsw.p_from_z(-zl,85.0*np.ones(zl.shape[0]))


def zmld_boyer(s, t, p):
    """
    Computes mixed layer depth, based on de Boyer Montégut et al., 2004.
    Parameters
    ----------
    s : array_like
        salinity [psu (PSS-78)]
    t : array_like
        temperature [℃ (ITS-90)]
    p : array_like
        pressure [db].
    Notes
    -----
    Based on density with fixed threshold criteria
    de Boyer Montégut et al., 2004. Mixed layer depth over the global ocean:
        An examination of profile data and a profile-based climatology.
        doi:10.1029/2004JC002378
    dataset for test and more explanation can be found at:
    http://www.ifremer.fr/cerweb/deboyer/mld/Surface_Mixed_Layer_Depth.php
    Codes based on : http://mixedlayer.ucsd.edu/
    """
    m = len(np.nonzero(~np.isnan(s))[0])

    if m <= 1:
        mldepthdens_mldindex = 0
        mldepthptemp_mldindex = 0
        return mldepthdens_mldindex, mldepthptemp_mldindex
    else:
        # starti = min(find((pres-10).^2==min((pres-10).^2)));
        starti = np.min(np.where((p - 10.0) ** 2 == np.min((p - 10.0) ** 2))[0])
        starti = 0
        pres = p[starti:m]
        sal = s[starti:m]
        temp = t[starti:m]

        pden = gsw.rho(sal, temp, 0) - 1000

        mldepthdens_mldindex = m - 1
        for i, pp in enumerate(pden):
            if np.abs(pden[starti] - pp) > 0.03:
                mldepthdens_mldindex = i
                break

        # Interpolate to exactly match the potential density threshold.
        presseg = [pres[mldepthdens_mldindex - 1], pres[mldepthdens_mldindex]]
        pdenseg = [
            pden[starti] - pden[mldepthdens_mldindex - 1],
            pden[starti] - pden[mldepthdens_mldindex],
        ]
        P = np.polyfit(presseg, pdenseg, 1)
        presinterp = np.linspace(presseg[0], presseg[1], 3)
        pdenthreshold = np.polyval(P, presinterp)

        # The potential density threshold MLD value:
        ix = np.max(np.where(np.abs(pdenthreshold) < 0.03)[0])
        mldepthdens_mldindex = presinterp[ix]

        # Search for the first level that exceeds the temperature threshold.
        mldepthptmp_mldindex = m - 1
        for i, tt in enumerate(temp):
            if np.abs(temp[starti] - tt) > 0.2:
                mldepthptmp_mldindex = i
                break

        # Interpolate to exactly match the temperature threshold.
        presseg = [pres[mldepthptmp_mldindex - 1], pres[mldepthptmp_mldindex]]
        tempseg = [
            temp[starti] - temp[mldepthptmp_mldindex - 1],
            temp[starti] - temp[mldepthptmp_mldindex],
        ]
        P = np.polyfit(presseg, tempseg, 1)
        presinterp = np.linspace(presseg[0], presseg[1], 3)
        tempthreshold = np.polyval(P, presinterp)

        # The temperature threshold MLD value:
        ix = np.max(np.where(np.abs(tempthreshold) < 0.2)[0])
        mldepthptemp_mldindex = presinterp[ix]

        return mldepthdens_mldindex, mldepthptemp_mldindex


