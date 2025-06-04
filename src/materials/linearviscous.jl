export liV!

@kernel inbounds = true function liV!(
    mp  ::DeviceParticle2D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT_1::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        nid  = attr.nid[ix]
        ηlin = attr.ηlin[nid]       # shear viscosity
        ηblk = attr.ηblk[nid]       # bulk viscosity

        # get strain rate for xx, yy, xy
        dϵxx = ΔT_1 * mp.ΔFs[ix, 1] 
        dϵyy = ΔT_1 * mp.ΔFs[ix, 4]
        dϵxy = ΔT_1 * (mp.ΔFs[ix, 2] + mp.ΔFs[ix, 3]) * T2(0.5)

        # pressure
        ϵvol = dϵxx + dϵyy
        p   = ϵvol/ηblk
        
        # deviatoric stress
        sxx = T2(2.0) * ηlin * dϵxx
        syy = T2(2.0) * ηlin * dϵyy
        sxy = T2(2.0) * ηlin * dϵxy
        
        # update stress
        mp.σij[ix, 1] = sxx + p
        mp.σij[ix, 2] = syy + p
        mp.σij[ix, 4] = sxy
        mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2]) * T2(0.5)
        mp.sij[ix, 1] = mp.σij[ix, 1] - mp.σm[ix]
        mp.sij[ix, 2] = mp.σij[ix, 2] - mp.σm[ix]
        mp.sij[ix, 4] = mp.σij[ix, 4]
    end
end

@kernel inbounds = true function liV!(
    mp  ::DeviceParticle3D{T1, T2},
    attr::  DeviceProperty{T1, T2},
    ΔT_1::T2
) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.np
        nid  = attr.nid[ix]
        ηlin = attr.ηlin[nid]
        
        # get strain rate for xx, yy, zz, xy, yz, zx
        dϵxx = ΔT_1 * mp.ΔFs[ix, 1] 
        dϵyy = ΔT_1 * mp.ΔFs[ix, 5]
        dϵzz = ΔT_1 * mp.ΔFs[ix, 9]
        dϵxy = ΔT_1 * (mp.ΔFs[ix, 2] + mp.ΔFs[ix, 4]) * T2(0.5)
        dϵyz = ΔT_1 * (mp.ΔFs[ix, 6] + mp.ΔFs[ix, 8]) * T2(0.5)
        dϵzx = ΔT_1 * (mp.ΔFs[ix, 3] + mp.ΔFs[ix, 7]) * T2(0.5)
        
        # pressure
        ϵvol = dϵxx + dϵyy
        p   = ϵvol/ηblk

        # deviatoric stress
        sxx = T2(2.0) * ηlin * dϵxx
        syy = T2(2.0) * ηlin * dϵyy
        szz = T2(2.0) * ηlin * dϵzz
        sxy = T2(2.0) * ηlin * dϵxy
        syz = T2(2.0) * ηlin * dϵyz
        szx = T2(2.0) * ηlin * dϵzx
        
        # update stress
        mp.σij[ix, 1] = sxx + p
        mp.σij[ix, 2] = syy + p
        mp.σij[ix, 3] = szz + p
        mp.σij[ix, 4] = sxy
        mp.σij[ix, 5] = syz
        mp.σij[ix, 6] = szx
        mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.sij[ix, 1] = mp.σij[ix, 1] - mp.σm[ix]
        mp.sij[ix, 2] = mp.σij[ix, 2] - mp.σm[ix]
        mp.sij[ix, 3] = mp.σij[ix, 3] - mp.σm[ix]
        mp.sij[ix, 4] = mp.σij[ix, 4]
        mp.sij[ix, 5] = mp.σij[ix, 5]
        mp.sij[ix, 6] = mp.σij[ix, 6]
    end
end
