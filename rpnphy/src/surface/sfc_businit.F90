!-------------------------------------- LICENCE BEGIN --------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ----------------------------

!/@*
subroutine sfc_businit(moyhr,ni,nk)
   use sfc_options
   use sfcbus_mod
   use svs_configs
   use class_configs
   implicit none
#include <arch_specific.hf>
   !@Object Establishes requirements in terms of variables in the 4 main buses
   !        (busent, busdyn, busper and busvol) for the TEB (Urban scheme).
   !@Arguments
   integer, intent(in) ::  moyhr,ni,nk !# horiz and vert dimensions
   !@Author M. Desgagne (Oct 1995)

   !@Revision
   ! 001      L. Spacek  (Aug 2010) - Complete rewrite
   ! 002      B. Dugas   (Oct 2010) - A few small corrections
   ! 003      L. Spacek  (Sep 2011) - Eliminate obsolete convection options
   ! 004      M. Abrahamowicz (May 2016) - Add SVS
   ! 005      K. Winger (ESCER/UQAM) (Feb 2017) - Add variables for lake fraction
   ! 005      K. Winger (ESCER/UQAM) (Jun 2020) - Variables for CLASS added
   !*@/

#include "phymkptr.hf"

   integer :: nmos
   integer :: alb_road, alb_roaden, alb_roof, alb_roofen, alb_wall, &
        alb_wallen, azim, bld, blden, bld_height, bld_heighten, can_hw_ratio, &
        d_road, d_roaden, d_roof, d_roofen, d_wall, d_wallen, emis_road, &
        emis_roaden, emis_roof, emis_roofen, emis_wall, emis_wallen, fvapliq, &
        fvapliqaf, g_road, &
        g_roof, g_town, g_wall, h_industry, h_industryen, h_road, h_roof, &
        h_town, h_traffic, h_trafficen, h_wall, hc_road, hc_roaden, hc_roof, &
        hc_roofen, hc_wall, hc_wallen, le_industry, le_industryen, le_road, &
        le_roof, le_town, le_traffic, le_trafficen, le_wall, nat, naten, pav, &
        paven, q_canyon, q_canyonen, rn_road, rn_roof, rn_town, rn_wall, &
        runofftot, runofftotaf, &
        sroad_alb, sroad_alben, sroad_emis, sroad_emisen, sroad_rho, &
        sroad_rhoen, sroad_scheme, sroad_t, sroad_ten, sroad_ts, sroad_tsen, &
        sroad_wsnow, sroad_wsnowen, sroof_alb, sroof_alben, sroof_emis, &
        sroof_emisen, sroof_rho, sroof_rhoen, sroof_scheme, sroof_t, &
        sroof_ten, sroof_ts, sroof_tsen, sroof_wsnow, sroof_wsnowen, &
        svf_road, svf_wall, t_canyon, t_canyonen, t_road, t_roaden, t_roof, &
        t_roofen, t_wall, t_wallen, tc_road, tc_roaden, tc_roof, tc_roofen, &
        tc_wall, tc_wallen, ti_bld, ti_blden, ti_road, ti_roaden, tsun, &
        u_canyon, wall_o_hor, wall_o_horen, ws_road, ws_roaden, ws_roof, &
        ws_roofen, &
        yradin, yradrfsun, yradrfshade, yutciin, yutcirfsun, &
        yutcirfshade, ytrfzt, ytrdzt, yurdzu, ywbgtrfsun, ywbgtrfshade, &
        yutcicin, yutcicsun, yutcicshade, yutcicrfsun, yutcicrfshade, &
        ytglbrfsun, ytglbrfshade, ytwetbrf, yq8, yq9, yq10, yq11, yq12, &
        yq13, &
        z0_road, z0_roaden, z0_roof, z0_roofen, z0_town, &
        z0_townen, zenith, emtw, alscatw, tsradtw
   integer :: acoef, alveg, bcoef, c1sat, c2ref, c3ref, clay, cveg, &
        eflux, emsvc, gamveg, husurf, hv, iceline, lai, melts,  &
        meltsr, pcoef, psn, psng, psnv, resa, rgl, rnet_s, rst, sand, &
        snoagen, snoalen, snoma, snoro, stomr, tsoil, vegf, &
        vegfrac, wfc, wsat, wsnow, wveg, wwilt
   integer :: cgsat, dsst, dtdiag, glacier, glsea0, &
        icedp, sfcwgt, skin_depth, skin_inc, snoal, snoden, &
        snodp, tglacier, tmice, tnolim, &
        twater, urban, &
        yradsun, yradshade, yutcisun, yutcishade, &
        ywbgtsun, ywbgtshade, ytglbsun, ytglbshade, ytwetb, yQ1, yQ2, &
        yq3, yq4, yq5, yq6, yq7, &
        z0en, z0veg, z0tveg, qdiagtyp, tdiagtyp, udiagtyp, vdiagtyp, &
        qdiagtypv, tdiagtypv, udiagtypv, vdiagtypv, tddiagtyp, tddiagtypv
   ! lake fields
   integer :: lakefr, frv_li, frv_lw, lakect, lakedepth, lakefice, lakehice, &
        lakehml, laketbot, laketice, laketmnw, laketp, laketransp, laketwml
   character(len=2) :: nm, nagg, nrow, ncg, ncv, ncvp, nicc, niccp
   character(len=3) :: ncvxcg, niccxcg
   !--------   FOR SVS -----------------
   character(len=2) :: ngl, nglp1, nstel, nstpl, iemib, iicel
   integer :: accevap,acroot, algr, alvl , alvh, avg_gwsol, clayen, co2i1, cvh, cvl, d50, d95, &
        deciduous, draindens, eg, emis, emisgr, emistg, emistgen, emisvh, emisvl, &
        er, etr, evergreen, &
        fbcof, frootd, gamvh, gamvl, grkef, hfluxsa, hfluxsv, &
        impervu, &
        khc, ksat, ksatc, laictem, laideci, laiva, laivf26, laivh, laivl, &
        latflaf, latflw, lesv, psi, psngrvl, psnvh, psnvha, &
        rcctem, resagr, resavg, resasa, resasv, resaef, rglvh, rglvl, &
        rnetsa, rnetsv, rsnowsa, &
        rsnowsv, rveg, sanden, skyview, slop, snodpl, snval, &
        snvden,  snvdp, snvma,  snvro, stomrvh, stomrvl, svs_wta, &
        tground, tsa, tsnavg, tsnow, tsnowveg, &
        tsvavg, tvege,  vegh, vegl, vegtrans, vgctem, &
        watflow, wsoilm, wfcdp, wsnv, &
        z0ha, z0mvh, z0mvl
   !--------   for SVS and CLASS ---------
   integer :: grksat, psisat, wfcint
   !--------   FOR CLASS -----------------
   integer :: alvs, alir, algwv, algwn, algdv, algdn, &
        ail, pai, alirc, alvsc, bbi, cdh, cdm, &
        cmai, delzw, evapo, fcanmx, fcovc, fcovcs, fcovg, fcovgs, firupaf, &
        flgg, flgs, flgv, fsgg, fsgs, fsgv, fsnow, fsolupaf, grkfac, &
        hcps, hevc, hevg, hevs, hfsc, hfsg, hfss, hmfc, hmfg, hmfn, &
        htc, htcc, htcs, huaircan, iveg, laimax, laimin, &
        mosfract, orgm, pcfc, pcpn, pclc, pcpg, psiga, &
        psigb, psiwlt, qa50, qfc, qfcf, qfcl, qfg, qfn, rib, &
        rofc, rofn, rovg, sdepth, &
        subflw, taircan, tbase, &
        tbasfl, tcs, thfc, thlmin, thlrat, thlret, thpor, tovrfl, &
        tpond, trunoff, tsno, tsubfl, tsurfsa, tveg, &
        veggro, vegma, vpda, vpdb, &
        wfsurf, wtrc, wtrg, wtrs, xdrain, xslope, zbotw, zoln, zpond, &
        zponden
   integer :: anis, are, excw, lbedr, leggw, slpgw, totw, wtnew

! Parameters for CLASS 3.6
   integer :: algwet, algdry
   !--------   FOR CTEM -----------------
   integer :: ailc, ailcb, ailcg, alirctm, allwacc, alswacc, alvsctm, &
        ancgvgac, ancsvgac, anndefct, annpcp, annsrpls, anpcpcur, anpecur, &
        anpotevp, aridity, bleafmas, bmasveg, burnvegf, cfluxcg, &
        cfluxcs, cmasvegc, co2conc, co2i1cg, co2i1cs, co2i2cg, co2i2cs, &
        colddayr, defctcur, defctmon, defmnr, dftcuryr, dryslen, dvdfcan, &
        extnprob, fcancmx, flhrloss, flinacc, flutacc, fsinacc, &
        fsnowacc, gavglai, gavgltms, gavgscms, gdd5, gdd5cur, geremort, &
        gleafmas, grwtheff, intrmort, lambda, lfstatur, &
        litrmass, lyrotmas, lystmmas, mlightng, nfcancmx, nppveg, &
        paic, pandayr, pfcancmx, pftexistr, pgleafmass, prbfrhuc, &
        preacc, pstemmass, rmatc, rmatctem, rmlcgvga, rmlcsvga, rootdpth, &
        rootmass, rothrlos, slai, slaic, soilcmas, &
        srpcuryr, srplscur, srplsmon, stemmass, stmhrlos, surmnr, &
        taaccgat, tbaraccgat, tbarcacc, tbarcsacc, tbargacc, tbargsacc, &
        tcanoaccgat, tcansacc, tcoldm, tcurm, thicecacc, thliqcacc, &
        thliqgacc, tmonthb, todfrac, twarmm, tymaxlai, uvaccgat, veghght, &
        vgbiomas, vvaccgat, wdmindex, zolnc
   integer :: afrleaf, afrroot, afrstem, autores, autresveg, burnarea, &
        colrate, dstcemls, dstcemls3, gpp, gppveg, grclarea, hetresveg, &
        hetrores, humiftrs, leaflitr, litrfall, litres, litresveg, lucemcom, &
        lucltrin, lucsocin, ltstatus, mortrate, nbp, nbpveg, nep, nepveg, npp, &
        probfire, rg, rgveg, rm, rml, rmlvegacc, rmr, rmrveg, rms, rmsveg, &
        roottemp, socres, socresveg, soilresp, tltrleaf, tltrroot, tltrstem, &
        vgbiomas_veg, wtstatus
   ! For CLASSIC
   integer :: grdhflx, soilcol, snowsize

 

   if (schmsol == 'SVS') then
      ! initialize levels for soil texture data
      call init_soil_text_levels()

      ! number of soil/"ground" layers in svs
      Write(ngl,'(i2)') nl_svs
      ! number of soil/"ground" layers PLUS 1
      Write(nglp1,'(i2)') nl_svs+1
      ! number of layer for ENTRY bus SVS clay and sand var.
      Write(nstel,'(i2)') nl_ste

      ! number of layer for PHYSICS bus SVS clay and sand var.
      Write(nstpl,'(i2)') nl_stp

   elseif (schmsol == 'CLASS') then

      ! number of special vegetation classes (CLASS)
      Write(ncv,'(i2)') CLASS_IC
      Write(ncvp,'(i2)') (CLASS_IC+1)

      ! number of ground layers (CLASS)
      Write(ncg,'(i2)') class_IG

      ! needed for rooting depth
      Write(ncvxcg,'(i3)') (CLASS_IC*class_IG)

     ! number of plant functional types (PFTs) in CTEM
     Write(nicc,'(i2)') CTEM_ICC

     ! additional parameters used when running CTEM
     Write(niccp,'(i2)') (CTEM_ICC+1)
     Write(ncvxcg,'(i3)') (CLASS_IC*class_IG)
     Write(niccxcg,'(i3)') (CTEM_ICC*class_IG)

   endif

   !---------------------------------------------------------------------
   write(nagg,'(i2)') nsurf+1

   !# nm is the number of mosaic 'layers'
   nmos = 0
   write(nm,'(i2)') nmos !#TODO: delete

   !# nl is the number of levels in sea ice
   write(nrow,'(i2)') nl

   !# EMIB is only a required input if isba_soil_emiss=='CLIMATO'
   iemib = '0'
   if (isba_soil_emiss=='CLIMATO') iemib = '1'

   !# ICEL is only a required input if icelac=.true.
   iicel = '0'
   if (icelac) iicel = '1'

   !#TODO: check if schmsol conditional
   PHYVAR2D1(cgsat,        'VN=cgsat        ;ON=6I  ;VD=thermal coef. at saturation                                          ;VB=p0')
   PHYVAR2D1(dsst,         'VN=dsst         ;ON=DSST;VD=warm layer diurnal SST increment                                     ;VB=p0')
   PHYVAR2D1(dtdiag,       'VN=dtdiag       ;ON=DLIM;VD=DeltaT at screen level of tdiaglim                                   ;VB=p0')
   PHYVAR2D1(glacier,      'VN=glacier      ;ON=2F  ;VD=continental ice fraction                                             ;VB=p1;IN=GA  ;MIN=0')
   PHYVAR2D1(glsea0,       'VN=glsea0       ;ON=GY  ;VD=sea ice fraction (unmodified)                                        ;VB=p1;IN=LG  ;MIN=0')
   PHYVAR2D1(icedp,        'VN=icedp        ;ON=I8  ;VD=sea ice thickness                                                    ;VB=p1        ;MIN=0')
   PHYVAR2D1(iceline,      'VN=iceline      ;ON=ICEL;VD=ice line                                                             ;VB=p'//iicel)
   PHYVAR2D1(lakefr,       'VN=lakefr       ;ON=FU  ;VD=lake fraction                                                        ;VB=p0        ;MIN=0')
   PHYVAR3D1(qdiagtyp,     "VN=qdiagtyp     ;ON=DQST;VD=screen level specific humidity for each sfc type   ; MIN=0 ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(qdiagtypv,    "VN=qdiagtypv    ;ON=DQSZ;VD=qdiagtyp for z0 vegetation-only ; MIN=0 ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR2D1(skin_depth,   'VN=skin_depth   ;ON=SDEP;VD=sea surface cold skin depth                                          ;VB=p0')
   PHYVAR2D1(skin_inc,     'VN=skin_inc     ;ON=SINC;VD=sea surface cold skin SST increment                                  ;VB=p0')
   PHYVAR3D1(snodp,        'VN=snodp        ;ON=SD  ;VD=snow depth                                     ;VS=A*'//nagg//'    ;VB=p1        ;MIN=0')
   PHYVAR3D1(tddiagtyp,    "VN=tddiagtyp    ;ON=TDST;VD=screen level dew temperature for each sfc type    ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tddiagtypv,   "VN=tddiagtypv   ;ON=TDSZ;VD=tddiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tdiagtyp,     "VN=tdiagtyp     ;ON=TT2M;VD=screen level temperature for each sfc type    ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tdiagtypv,    "VN=tdiagtypv    ;ON=TJSZ;VD=tdiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(tglacier,     'VN=tglacier     ;ON=I9  ;VD=glaciers temperature                           ;VS=A*2             ;VB=p1')
   PHYVAR3D1(tmice,        'VN=tmice        ;ON=I7  ;VD=sea ice temperature                            ;VS=A*'//nrow//'    ;VB=p1')
   PHYVAR2D1(tnolim,       'VN=tnolim       ;ON=TNOL;VD=screen level temp without max on gradient                            ;VB=p0')
   PHYVAR2D1(twater,       'VN=twater       ;ON=TM  ;VD=sea surface temperature                                              ;VB=p1')
   PHYVAR3D1(udiagtyp,     "VN=udiagtyp     ;ON=UDST;VD=screen level U-wind for each sfc type         ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(udiagtypv,    "VN=udiagtypv    ;ON=UDSZ;VD=udiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vdiagtyp,     "VN=vdiagtyp     ;ON=VDST;VD=screen level V-wind for each sfc type         ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vdiagtypv,    "VN=vdiagtypv    ;ON=VDSZ;VD=vdiagtyp for z0 vegetation-only ; VS=A*"//nagg//" ; VB=v0")
   PHYVAR3D1(vegf,         'VN=vegf         ;ON=2V  ;VD=vegetation fractions                           ;VS=A*26            ;VB=p1        ;MIN=0; IN=VF')
   PHYVAR3D1(vegfrac,      'VN=vegfrac      ;ON=K1  ;VD=vegetation fraction                            ;VS=A*5             ;VB=p0        ;MIN=0')
   PHYVAR2D1(urban,        'VN=urban        ;ON=UR  ;VD=urban mask                                                           ;VB=p0')
 
   PHYVAR3DC(yradsun,      "VN=yradsun      ;ON=RTSU;VD=MRT in the exposed sunny street (K)                       ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yradshade,    "VN=yradshade    ;ON=RTHD;VD=MRT in the shaded street (K)                              ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yutcisun,     "VN=yutcisun     ;ON=DXSU;VD=UTCI in the exposed sunny street (C)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yutcishade,   "VN=yutcishade   ;ON=DXHD;VD= UTCI in the shaded street (C)                            ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ywbgtsun,     "VN=ywbgtsun     ;ON=GXSU;VD= WBGT in the exposed sunny street (C)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ywbgtshade,   "VN=ywbgtshade   ;ON=GXHD;VD= WBGT in the shaded street (C)                            ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytglbsun,     "VN=ytglbsun     ;ON=GTSU;VD=TGlobe in the exposed sunny street (K)                    ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytglbshade,   "VN=ytglbshade   ;ON=GTHD;VD=TGlobe in the shaded street (K)                           ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(ytwetb,       "VN=ytwetb       ;ON=WBT ;VD=Wet-Bulb Temp at zt above the ground (K)                  ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ1,          "VN=yQ1          ;ON=QSSU;VD= Contribution of direct solar rad (W/m2)                  ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ2,          "VN=yQ2          ;ON=QSSK;VD= Contribution of sky SW rad (W/m2)                        ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ3,          "VN=yQ3          ;ON=QLSK;VD= Contribution of sky LW rad (W/m2)                        ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ4,          "VN=yQ4          ;ON=QSRD;VD= Contribution of ground SW rad (W/m2)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ5,          "VN=yQ5          ;ON=QLRD;VD= Contribution of ground LW rad (W/m2)                     ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ6,          "VN=yQ6          ;ON=QSWL;VD= Contribution of facet SW rad (W/m2)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)
   PHYVAR3DC(yQ7,          "VN=yQ7          ;ON=QLWL;VD= Contribution of facet LW rad (W/m2)                      ; VS=A*"//nagg//" ; VB=v0", thermal_stress)

   PHYVAR2D1(z0veg,        'VN=z0veg        ;ON=ZVG2;VD=vegetation roughness length                                          ;VB=p1;MIN=0')
   PHYVAR2D1(z0tveg,       'VN=z0tveg       ;ON=ZVGT;VD=thermodynamic vegetation roughness length                            ;VB=p0;MIN=0')
   if (z0veg_only) then
      PHYVAR2D1(z0en,         'VN=z0en         ;ON=2B  ;VD=roughness length (E; vegetation only)                             ;VB=e1; IN=ZVG2; MIN=0')
   else
      PHYVAR2D1(z0en,         'VN=z0en         ;ON=2B  ;VD=roughness length (E)                                              ;VB=e1; IN=ZP')
   endif

   if (moyhr > 0) &
        PHYVAR2D1(insmavg,   'VN=insmavg      ;ON=IMAV;VD=integrated soil moist avg over last moyhr hrs                      ;VB=p0        ;MIN=0')

   ! Common surface fields
   PHYVAR2D1(fvapliq,      'VN=fvapliq      ;ON=HFLQ;VD=surf. evaporation (kg/m2 or mm)                                   ;VB=p0')
   PHYVAR2D1(fvapliqaf,    'VN=fvapliqaf    ;ON=AHFL;VD=accum. surf. evaporation (HFLQ) (kg/m2 or mm)                     ;VB=p0')
   PHYVAR3D1(runofftot,    'VN=runofftot    ;ON=TRUN;VD=total surface runoff                           ;VS=A*'//nagg//' ;VB=v0')
   PHYVAR3D1(runofftotaf,  'VN=runofftotaf  ;ON=TRAF;VD=accum. of total surface runoff                 ;VS=A*'//nagg//' ;VB=p0')
   PHYVAR2D1(snoma,        'VN=snoma        ;ON=I5  ;VD=snow mass                                                         ;VB=p0        ;MIN=0')

   IF_ISBA: if (schmsol == 'ISBA') then
      PHYVAR2D1(acoef,        'VN=acoef        ;ON=1I  ;VD=a coef. in wgeq                                                   ;VB=p0')
      PHYVAR2D1(alveg,        'VN=alveg        ;ON=AX  ;VD=visible canopy albedo                                             ;VB=p0        ;MIN=0')
      PHYVAR2D1(bcoef,        'VN=bcoef        ;ON=1G  ;VD=slope of retention curve                                          ;VB=p0')
      PHYVAR2D1(c1sat,        'VN=c1sat        ;ON=3I  ;VD=c1 coef. at saturation                                            ;VB=p0')
      PHYVAR2D1(c2ref,        'VN=c2ref        ;ON=4I  ;VD=reference value of c2                                             ;VB=p0')
      PHYVAR2D1(c3ref,        'VN=c3ref        ;ON=5I  ;VD=drainage coef. to deeper soil                                     ;VB=p0')
      PHYVAR3D1(clay,         'VN=clay         ;ON=J2  ;VD=percentage of clay in soil                     ;VS=A*3            ;VB=p1        ;MIN=0')
      PHYVAR2D1(cveg,         'VN=cveg         ;ON=CV  ;VD=thermal coefficient for canopy                                    ;VB=p0')
      PHYVAR2D1(drain,        'VN=drain        ;ON=DR  ;VD=water drainage at bottom of soil layer                            ;VB=p0')
      PHYVAR2D1(drainaf,      'VN=drainaf      ;ON=O1  ;VD=accum. of base drainage                                           ;VB=p0')
      PHYVAR2D1(eflux,        'VN=eflux        ;ON=4F  ;VD=specific hum. flux (=-alfaq)                                      ;VB=v0')
      PHYVAR2D1(emsvc,        'VN=emsvc        ;ON=EMIB;VD=broadband climatological ground-veg emissivity isba               ;VB=p'//iemib//';MIN=0')
      PHYVAR2D1(gamveg,       'VN=gamveg       ;ON=GG  ;VD=stomatal resistance parameter                                     ;VB=p0')
      PHYVAR2D1(husurf,       'VN=husurf       ;ON=FH  ;VD=spec. humid. of the surface                                       ;VB=v0        ;MIN=0')
      PHYVAR2D1(hv,           'VN=hv           ;ON=HV  ;VD=relative humidity of veg. canopy                                  ;VB=v0        ;MIN=0')
      PHYVAR2D1(isoil,        'VN=isoil        ;ON=I2  ;VD=soil volumetric ice contents                                      ;VB=p1        ;MIN=0')
      PHYVAR2D1(lai,          'VN=lai          ;ON=J4  ;VD=leaf area index                                                   ;VB=p0        ;MIN=0')
      PHYVAR2D1(leg,          'VN=leg          ;ON=L2  ;VD=latent heat flux over bare grnd                                   ;VB=v0')
      PHYVAR2D1(legaf,        'VN=legaf        ;ON=O5  ;VD=accum. of bare ground LE flux                                     ;VB=p0')
      PHYVAR2D1(ler,          'VN=ler          ;ON=LR  ;VD=latent heat flux from leaves                                      ;VB=v0')
      PHYVAR2D1(leraf,        'VN=leraf        ;ON=O6  ;VD=accum. of direct veg LE flux                                      ;VB=p0')
      PHYVAR2D1(les,          'VN=les          ;ON=LS  ;VD=latent heat flux over snow                                        ;VB=v0')
      PHYVAR2D1(lesaf,        'VN=lesaf        ;ON=O7  ;VD=accum. of sublimation from snow                                   ;VB=p0')
      PHYVAR2D1(letr,         'VN=letr         ;ON=LT  ;VD=latent heat of evapotransp.                                       ;VB=v0')
      PHYVAR2D1(letraf,       'VN=letraf       ;ON=O8  ;VD=accum. of veg. transpiration                                      ;VB=p0')
      PHYVAR2D1(lev,          'VN=lev          ;ON=LV  ;VD=latent heat flux over vegetation                                  ;VB=v0')
      PHYVAR2D1(levaf,        'VN=levaf        ;ON=O9  ;VD=accum. of evaporation from veg.                                   ;VB=p0')
      PHYVAR2D1(melts,        'VN=melts        ;ON=MLTS;VD=accum. snow melting (kg/m2)                                       ;VB=p0')
      PHYVAR2D1(meltsr,       'VN=meltsr       ;ON=MLTR;VD=accum. snow melting due to rain (kg/m2)                           ;VB=p0')
      PHYVAR2D1(overfl,       'VN=overfl       ;ON=RO  ;VD=overland runoff                                                   ;VB=v0')
      PHYVAR2D1(overflaf,     'VN=overflaf     ;ON=N0  ;VD=accum. of surface runoff                                          ;VB=p0')
      PHYVAR2D1(pcoef,        'VN=pcoef        ;ON=7I  ;VD=p coef. in wgeq                                                   ;VB=p0')
      PHYVAR2D1(psn,          'VN=psn          ;ON=5P  ;VD=fraction of the grid covered by snow                              ;VB=v0        ;MIN=0')
      PHYVAR2D1(psng,         'VN=psng         ;ON=3P  ;VD=fraction of bare ground covered by snow                           ;VB=v0        ;MIN=0')
      PHYVAR2D1(psnv,         'VN=psnv         ;ON=4P  ;VD=fraction of vegetation covered by snow                            ;VB=v0        ;MIN=0')
      PHYVAR2D1(resa,         'VN=resa         ;ON=RD  ;VD=aerodynamic resistance                                            ;VB=p0')
      PHYVAR2D1(rgl,          'VN=rgl          ;ON=RG  ;VD=parameter stomatal resistance                                     ;VB=p0')
      PHYVAR2D1(rnet_s,       'VN=rnet_s       ;ON=NR  ;VD=net radiation (soil only)                                         ;VB=v0')
      PHYVAR2D1(rootdp,       'VN=rootdp       ;ON=D2  ;VD=rooting soil depth                                                ;VB=p0')
      PHYVAR2D1(rst,          'VN=rst          ;ON=R1  ;VD=stomatal resistance                                               ;VB=v0')
      PHYVAR3D1(sand,         'VN=sand         ;ON=J1  ;VD=percentage of sand in soil                     ;VS=A*3          ;VB=p1        ;MIN=0')
      if (snoalb_anl) then
         PHYVAR2D1(snoalen,      'VN=snoalen      ;ON=5H  ;VD=snow albedo (E)                                                ;VB=e1;IN=I6  ;MIN=0')
      else
         PHYVAR2D1(snoagen,      'VN=snoagen      ;ON=3H  ;VD=age of snow (E)                                                ;VB=e1;IN=XA  ;MIN=0')
      endif
      PHYVAR3D1(snoal,        'VN=snoal        ;ON=I6  ;VD=albedo of snow                                 ;VS=A@'//nm//'   ;VB=p0        ;MIN=0')
      PHYVAR3D1(snoden,       'VN=snoden       ;ON=DN  ;VD=snow density in kg/m3                          ;VS=A@'//nm//'   ;VB=p0        ;MIN=0; IN=DN0')
      PHYVAR2D1(snoro,        'VN=snoro        ;ON=7S  ;VD=relative snow density                                             ;VB=p1;IN=DN  ;MIN=0')
      PHYVAR2D1(stomr,        'VN=stomr        ;ON=RS  ;VD=minimum stomatal resistance                                       ;VB=p0')
      PHYVAR3D1(tsoil,        'VN=tsoil        ;ON=I0  ;VD=surface and soil temperatures                  ;VS=A*2          ;VB=p1')
      PHYVAR2D1(wfc,          'VN=wfc          ;ON=J5  ;VD=vol. water content at field cap.                                  ;VB=p0        ;MIN=0')
      PHYVAR2D1(wflux,        'VN=wflux        ;ON=M8  ;VD=water flux from surface to atm.                                   ;VB=v0')
      PHYVAR2D1(wfluxaf,      'VN=wfluxaf      ;ON=N7  ;VD=acc. of soil surface upward water flux                            ;VB=p0')
      PHYVAR2D1(wsat,         'VN=wsat         ;ON=J6  ;VD=vol. water content at saturation                                  ;VB=p0        ;MIN=0')
      PHYVAR3D1(wsnow,        'VN=wsnow        ;ON=I4  ;VD=water in the snow pack                         ;VS=A@'//nm//'   ;VB=p1        ;MIN=0')
      PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=I1  ;VD=soil volumetric water contents                 ;VS=A*2          ;VB=p1        ;MIN=0')
      PHYVAR2D1(wveg,         'VN=wveg         ;ON=I3  ;VD=water retained on the vegetation                                  ;VB=p1        ;MIN=0')
      PHYVAR2D1(wwilt,        'VN=wwilt        ;ON=J7  ;VD=vol. water cont. at wilting pt.                                   ;VB=p0        ;MIN=0')
   endif IF_ISBA

   IF_SVS: if (schmsol == 'SVS') then
! check/add min values !!!
      PHYVAR2D1(accevap,      'VN=accevap      ;ON=ACWF;VD=accum. of actual surf. evap. (kg/m2 or mm)                        ;VB=p0')
      PHYVAR3D1(acroot,       'VN=acroot       ;ON=ACRT;VD=active fraction of roots in soil layer         ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(algr,         'VN=algr         ;ON=ALGR;VD=visible albedo for bare ground                                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(alvh,         'VN=alvh         ;ON=ALVH;VD=visible canopy albedo for high vegetation only                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(alvl,         'VN=alvl         ;ON=ALVL;VD=visible canopy albedo for low vegetation only                     ;VB=p0        ;MIN=0')
      PHYVAR2D1(avg_gwsol,    'VN=avg_gwsol    ;ON=AGWS;VD=average soil moisture stress term                                 ;VB=p0        ;MIN=0')
      PHYVAR3D1(bcoef,        'VN=bcoef        ;ON=1G  ;VD=slope of retention curve                       ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(clay,         'VN=clay         ;ON=J2  ;VD=percentage of clay in soil                     ;VS=A*'//nstpl//';VB=p0        ;MIN=0')
      PHYVAR3D1(clayen,       'VN=clayen       ;ON=2H  ;VD=perc. of clay in soil (E)                      ;VS=A*'//nstel//';VB=e1;IN=J2  ;MIN=0')
      PHYVAR3D1(co2i1,        'VN=co2i1        ;ON=CO3 ;VD=CO2 CONCENTRATION   CTEM                       ;VS=A*9          ;VB=p0')
      PHYVAR2D1(cveg,         'VN=cveg         ;ON=CV  ;VD=thermal coefficient for canopy                                    ;VB=p0')
      PHYVAR2D1(cvh,          'VN=cvh          ;ON=CVH ;VD=thermal coefficient for canopy of high veg                        ;VB=p0')
      PHYVAR2D1(cvl,          'VN=cvl          ;ON=CVL ;VD=thermal coefficient for canopy of low veg                         ;VB=p0')
      PHYVAR2D1(d50,          'VN=d50          ;ON=d50 ;VD=depth[m] above which 50% of roots are located                     ;VB=p0')
      PHYVAR2D1(d95,          'VN=d95          ;ON=d95 ;VD=depth[m] above which 95% of roots are located                     ;VB=p0')
      PHYVAR2D1(deciduous,    'VN=deciduous    ;ON=DECI;VD=frac. of high veg. that is deciduous.                             ;VB=p0')
      PHYVAR3D1(drain,        'VN=drain        ;ON=DR  ;VD=water drainage in deep soil layers             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(drainaf,      'VN=drainaf      ;ON=O1  ;VD=accum. of base drainage                                           ;VB=p0')
      PHYVAR2D1(draindens,    'VN=draindens    ;ON=DRND;VD=drainage density (m/m2)                                           ;VB=p1')
      PHYVAR2D1(eflux,        'VN=eflux        ;ON=EFLX;VD=specific hum. flux (=-alfaq)                                      ;VB=v0')
      PHYVAR2D1(eg,           'VN=eg           ;ON=EG  ;VD=evapo. rate over bare grnd(no frac)                               ;VB=v0')
      PHYVAR2D1(emis,         'VN=emis         ;ON=EMI1;VD=emissivity of nat surface                                         ;VB=p0')
      PHYVAR2D1(emisgr,       'VN=emisgr       ;ON=EMGR;VD=emissivity of bare ground                                         ;VB=p0')
      PHYVAR2D1(emistg,       'VN=emistg       ;ON=EMTG;VD=emissivity land surface with no snow (read-in)                    ;VB=p0')
      if (read_emis) &
           PHYVAR2D1(emistgen,     'VN=emistgen     ;ON=ETG1;VD=avg. emissivity land surface with no snow (E)       ;VB=e1;IN=EMIB;MIN=0')
      PHYVAR2D1(emisvh,       'VN=emisvh       ;ON=EMVH;VD=emissivity of high vegetation                                     ;VB=p0')
      PHYVAR2D1(emisvl,       'VN=emisvl       ;ON=EMVL;VD=emissivity of low vegetation                                      ;VB=p0')
      PHYVAR2D1(er,           'VN=er           ;ON=ER  ;VD=evapo rate from leaves(no frac)                                   ;VB=v0')
      PHYVAR2D1(etr,          'VN=etr          ;ON=ETR ;VD=evapotranspiration rate (no frac)                                 ;VB=v0')
      PHYVAR2D1(evergreen,    'VN=evergreen    ;ON=EVER;VD=frac. of high veg. that is evergreen                              ;VB=p0')
      PHYVAR3D1(fbcof,        'VN=fbcof        ;ON=3G  ;VD=parameter derived from bcoef                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(frootd,       'VN=frootd       ;ON=FRTD;VD=deep soil layer root density                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(gamvh,        'VN=gamvh        ;ON=GGVH;VD=stomatal resistance parameter for high veg                        ;VB=p0')
      PHYVAR2D1(gamvl,        'VN=gamvl        ;ON=GGVL;VD=stomatal resistance parameter for low veg                         ;VB=p0')
      PHYVAR2D1(grkef,        'VN=grkef        ;ON=GKE; VD=WATDR parameter                                                   ;VB=p0')
      PHYVAR3D1(grksat,       'VN=grksat       ;ON=GKS  ;VD=sat. horiz. soil hydraulic conductivity       ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(hfluxsa,      'VN=hfluxsa      ;ON=HFSA;VD=sensible heat flux (snow only)                                    ;VB=p0')
      PHYVAR2D1(hfluxsv,      'VN=hfluxsv      ;ON=HFSV;VD=sensible heat flux (snow under veg. only)                         ;VB=p0')
      PHYVAR2D1(husurf,       'VN=husurf       ;ON=FH  ;VD=spec. humid. of the surface                                       ;VB=v0        ;MIN=0')
      PHYVAR2D1(hv,           'VN=hv           ;ON=HV  ;VD=relative humidity of veg. canopy                                  ;VB=v0        ;MIN=0')
      PHYVAR2D1(impervu,      'VN=impervu      ;ON=IMPU;VD=frac. of land sfc considered impervious (urban)                   ;VB=p0')
      PHYVAR3D1(isoil,        'VN=isoil        ;ON=ISOL;VD=soil volumetric ice contents per layer         ;VS=A*'//ngl//'  ;VB=p1        ;MIN=0')
      PHYVAR3D1(khc,          'VN=khc          ;ON=KHC ;VD=soil hydraulic conductivity                    ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(ksat,         'VN=ksat         ;ON=KSAT  ;VD=sat. soil hydraulic conductivity             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(ksatc,        'VN=ksatc        ;ON=KSTC  ;VD=corrected sat. soil hydraulic conductivity   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(laictem,      'VN=laictem      ;ON=LC    ;VD=vegetation LAI for 9 CTEM plant classes      ;VS=A*9          ;VB=p0')
      PHYVAR2D1(laideci,      'VN=laideci      ;ON=LAID;VD=leaf area index for high deciduous veg. only                      ;VB=p0')
      PHYVAR2D1(laiva,        'VN=laiva        ;ON=LAIA;VD=avg. leaf area index seen from atm.                               ;VB=p0')
      PHYVAR3D1(laivf26,      'VN=laivf26      ;ON=LAVF;VD=lai for each vf class times fractipm           ;VS=A*26         ;VB=p0')
      PHYVAR2D1(laivh,        'VN=laivh        ;ON=LAIH;VD=leaf area index for high vegetation only                          ;VB=p0')
      PHYVAR2D1(laivl,        'VN=laivl        ;ON=LAIL;VD=leaf area index for low vegetation only                           ;VB=p0')
      PHYVAR2D1(latflaf,      'VN=latflaf      ;ON=ALAT;VD=Accum. of LATF at all levels (kg/m2 = mm)                         ;VB=p0')
      PHYVAR3D1(latflw,       'VN=latflw       ;ON=LATF;VD=Lateral flow                                   ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(leg,          'VN=leg          ;ON=L2  ;VD=latent heat flux over bare grnd                                   ;VB=v0')
      PHYVAR2D1(ler,          'VN=ler          ;ON=LR  ;VD=latent heat flux from leaves                                      ;VB=v0')
      PHYVAR2D1(les,          'VN=les          ;ON=LS  ;VD=latent heat flux over snow                                        ;VB=v0')
      PHYVAR2D1(lesv,         'VN=lesv         ;ON=LSV ;VD=latent heat flux over snow-under-veg                              ;VB=v0')
      PHYVAR2D1(letr,         'VN=letr         ;ON=LT  ;VD=latent heat of evapotransp.                                       ;VB=v0')
      PHYVAR2D1(lev,          'VN=lev          ;ON=LV  ;VD=latent heat flux over vegetation                                  ;VB=v0')
      PHYVAR2D1(melts,        'VN=melts        ;ON=MLTS;VD=accum. snow melting (kg/m2)                                       ;VB=p0')
      PHYVAR2D1(meltsr,       'VN=meltsr       ;ON=MLTR;VD=accum. snow melting due to rain (kg/m2)                           ;VB=p0')
      PHYVAR3D1(psi,          'VN=psi          ;ON=PSI ;VD=soil water suction                             ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR3D1(psisat,       'VN=psisat       ;ON=D5  ;VD=sat. soil water suction                        ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(psngrvl,      'VN=psngrvl      ;ON=PSGL;VD=frac. of bare soil &/or low veg. cov. by snow                     ;VB=v0')
      PHYVAR2D1(psnvh,        'VN=psnvh        ;ON=PSVH;VD=fraction of high vegetation covered by snow                       ;VB=p0')
      PHYVAR2D1(psnvha,       'VN=psnvha       ;ON=PSVA;VD=frac. of high veg. covered by snow from atm.                      ;VB=p0')
      PHYVAR2D1(rcctem,       'VN=rcctem       ;ON=RCC ;VD=stomatal resistance CTEM                                          ;VB=p0')
      PHYVAR2D1(resagr,       'VN=resagr       ;ON=RSGR;VD=aerodynamic resistance over bare ground                           ;VB=p0')
      PHYVAR2D1(resavg,       'VN=resavg       ;ON=RSVG;VD=aerodynamic resistance over veget.                                ;VB=p0')
      PHYVAR2D1(resasa,       'VN=resasa       ;ON=RSSA;VD=aerodynamic resistance for snow on bg/low veg                     ;VB=p0')
      PHYVAR2D1(resasv,       'VN=resasv       ;ON=RSSV;VD=aerodynamic resistance for snow under high veg                    ;VB=p0')
      PHYVAR2D1(resaef,       'VN=resaef       ;ON=RSEF;VD=effective aerodynamic resistance for SVS land tile                ;VB=p0')
      PHYVAR2D1(rglvh,        'VN=rglvh        ;ON=RGVH;VD=parameter stomatal resistance for high veg                        ;VB=p0')
      PHYVAR2D1(rglvl,        'VN=rglvl        ;ON=RGVL;VD=parameter stomatal resistance for low veg                         ;VB=p0')
      PHYVAR2D1(rnet_s,       'VN=rnet_s       ;ON=NR  ;VD=net radiation (soil only)                                         ;VB=v0')
      PHYVAR2D1(rnetsa,       'VN=rnetsa       ;ON=RNSA;VD=net radiation (snow only)                                         ;VB=p0')
      PHYVAR2D1(rnetsv,       'VN=rnetsv       ;ON=RNSV;VD=net radiation (snow under veg. only)                              ;VB=p0')
      PHYVAR2D1(rootdp,       'VN=rootdp       ;ON=D2  ;VD=rooting soil depth                                                ;VB=p0')
      PHYVAR2D1(rsnowsa,      'VN=rsnowsa      ;ON=RSA ;VD=liquid water out of the snow pack                                 ;VB=p0')
      PHYVAR2D1(rsnowsv,      'VN=rsnowsv      ;ON=RSV ;VD=liquid water out of the snow-under-veg pack                       ;VB=p0')
      PHYVAR2D1(rst,          'VN=rst          ;ON=R1  ;VD=stomatal resistance                                               ;VB=v0')
      PHYVAR2D1(rveg,         'VN=rveg         ;ON=RVG ;VD=runoff from the vegetation (mm/s)                                 ;VB=p0')
      PHYVAR3D1(sand,         'VN=sand         ;ON=J1  ;VD=percentage of sand in soil                     ;VS=A*'//nstpl//';VB=p0')
      PHYVAR3D1(sanden,       'VN=sanden       ;ON=2G  ;VD=perc. of sand in soil (E)                      ;VS=A*'//nstel//';VB=e1;IN=J1  ;')
      PHYVAR2D1(skyview,      'VN=skyview      ;ON=SVF ;VD=sky view factor for tall vegetation                               ;VB=p0')
      PHYVAR2D1(slop,         'VN=slop         ;ON=SLOP;VD=average maximum subgrid-scale topo slope (nil)                    ;VB=p1')
      PHYVAR2D1(snoal,        'VN=snoal        ;ON=SNAL;VD=snow-over-low-veg/bare-ground albedo                              ;VB=p1        ;MIN=0')
      PHYVAR2D1(snoden,       'VN=snoden       ;ON=SNDN;VD=snow-over-low-veg/bare-ground density in kg/m3                    ;VB=p1        ;MIN=0')
      PHYVAR2D1(snodpl,       'VN=snodpl       ;ON=SNDP;VD=snow-over-low-veg/bare-ground depth                               ;VB=p1        ;MIN=0')
      PHYVAR2D1(snoro,        'VN=snoro        ;ON=SNDR;VD=snow-over-low-veg/bare-ground relative density                    ;VB=p0        ;MIN=0')
      PHYVAR2D1(snval,        'VN=snval        ;ON=SVAL;VD=snow-under-high-veg albedo                                        ;VB=p1')
      PHYVAR2D1(snvden,       'VN=snvden       ;ON=SVDN;VD=snow-under-high-veg density in kg/m3                              ;VB=p1')
      PHYVAR2D1(snvdp,        'VN=snvdp        ;ON=SVDP;VD=snow-under-high-veg depth                                         ;VB=p1')
      PHYVAR2D1(snvma,        'VN=snvma        ;ON=SVM ;VD=snow-under-high-veg mass                                          ;VB=p0')
      PHYVAR2D1(snvro,        'VN=snvro        ;ON=SVDR;VD=snow-under-high-veg relative density                              ;VB=p0')
      PHYVAR2D1(stomrvh,      'VN=stomrvh      ;ON=RSVH;VD=min. stomatal resistance for high vegetation                      ;VB=p0')
      PHYVAR2D1(stomrvl,      'VN=stomrvl      ;ON=RSVL;VD=min. stomatal resistance for low vegetation                       ;VB=p0')
      PHYVAR3D1(svs_wta,      'VN=svs_wta      ;ON=SVSW;VD=weight for svs used in aggregation **FROM SPACE;VS=A*5          ;VB=p0')
      PHYVAR3D1(tground,      'VN=tground      ;ON=TGR ;VD=skin and mean ground temp.                     ;VS=A*2          ;VB=p1')
      PHYVAR2D1(tsa,          'VN=tsa          ;ON=TSA ;VD=skin temp. of land surface as seen from atm                       ;VB=p0')
      PHYVAR2D1(tsnavg,       'VN=tsnavg       ;ON=ATSN;VD=snow-low-veg/bare-grnd avg temp. for melt/freez                   ;VB=p0')
      PHYVAR3D1(tsnow,        'VN=tsnow        ;ON=TSN ;VD=snow-low-veg/bare-grnd skin and mean temp.     ;VS=A*2          ;VB=p1')
      PHYVAR3D1(tsnowveg,     'VN=tsnowveg     ;ON=TSNV;VD=snow-under-high-veg skin and mean temp.        ;VS=A*2          ;VB=p1')
      PHYVAR2D1(tsvavg,       'VN=tsvavg       ;ON=ATSV;VD=snow-under-high-veg avg temp. for melt/freez                      ;VB=p0')
      PHYVAR3D1(tvege,        'VN=tvege        ;ON=TVG ;VD=skin and mean vegetation temp.                 ;VS=A*2          ;VB=p1')
      PHYVAR2D1(vegh,         'VN=vegh         ;ON=VEGH;VD=fraction of grid covered by high vegetation                       ;VB=p0')
      PHYVAR2D1(vegl,         'VN=vegl         ;ON=VEGL;VD=fraction of grid covered by low vegetation                        ;VB=p0')
      PHYVAR2D1(vegtrans,     'VN=vegtrans     ;ON=VGTR;VD=transmissivity of tall vegetation                                 ;VB=p0')
      PHYVAR3D1(vgctem,       'VN=vgctem       ;ON=VGCT;VD=CTEM vegetation type fractions                 ;VS=A*9          ;VB=p0')
      PHYVAR3D1(watflow,      'VN=watflow      ;ON=WFL ;VD=waterflow between layers                       ;VS=A*'//nglp1//';VB=p0')
      PHYVAR3D1(wfc,          'VN=wfc          ;ON=WFC ;VD=vol. water content at field cap.               ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wfcdp,        'VN=wfcdp        ;ON=WFCD;VD=vol. water content at field cap. at lowst layer                   ;VB=p0')
      PHYVAR3D1(wfcint,       'VN=wfcint       ;ON=WFCI;VD=water content at field capacity along slope    ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wflux,        'VN=wflux        ;ON=M8  ;VD=water flux from surface to atm.                                   ;VB=v0')
      PHYVAR2D1(wfluxaf,      'VN=wfluxaf      ;ON=N7  ;VD=acc. of soil surface upward water flux                            ;VB=p0')
      PHYVAR3D1(wsat,         'VN=wsat         ;ON=WSAT;VD=vol. water content at saturation               ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(wsnow,        'VN=wsnow        ;ON=WSN ;VD=water in low-veg/bare-grnd snowpack                               ;VB=p1        ;MIN=0')
      PHYVAR2D1(wsnv,         'VN=wsnv         ;ON=WSV ;VD=water in under-high-veg snowpack                                  ;VB=p1        ;MIN=0')
      PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=WSOL;VD=soil volm water content per layer              ;VS=A*'//ngl//'  ;VB=p1')
      PHYVAR2D1(wsoilm,       'VN=wsoilm       ;ON=WSLM;VD=mean soil volm watr cont for the whole column                     ;VB=p0')
      PHYVAR2D1(wveg,         'VN=wveg         ;ON=WVEG;VD=water retained on the vegetation                                  ;VB=p1        ;MIN=0')
      PHYVAR3D1(wwilt,        'VN=wwilt        ;ON=WWLT;VD=vol. water cont. at wilting pt.                ;VS=A*'//ngl//'  ;VB=p0')
      PHYVAR2D1(z0ha,         'VN=z0ha         ;ON=Z0HA;VD=thermal roughness for snowless veg.                               ;VB=p0')
      PHYVAR2D1(z0mvh,        'VN=z0mvh        ;ON=Z0VH;VD=local mom roughness length for high veg.                          ;VB=p0')
      PHYVAR2D1(z0mvl,        'VN=z0mvl        ;ON=Z0VL;VD=local mom roughness length for low veg.                           ;VB=p0')
   endif IF_SVS


   IF_CLASS: if (schmsol == 'CLASS') then
      PHYVAR2D1(alvs,         'VN=alvs         ;ON=ALVS;VD=Total visible albedo of land surface                              ;VB=v0')
      PHYVAR2D1(alir,         'VN=alir         ;ON=ALIR;VD=Total near-infrared albedo of land surface                        ;VB=v0')
      PHYVAR2D1(algdn,        'VN=algdn        ;ON=ALDN;VD=near-ir albedo of dry soil                                        ;VB=p0')
      PHYVAR2D1(algdv,        'VN=algdv        ;ON=ALDV;VD=visible albedo of dry soil                                        ;VB=p0')
      PHYVAR2D1(algwn,        'VN=algwn        ;ON=ALWN;VD=near-ir albedo of wet soil                                        ;VB=p0')
      PHYVAR2D1(algwv,        'VN=algwv        ;ON=ALWV;VD=visible albedo of wet soil                                        ;VB=p0')
      PHYVAR3D1(alirc,        'VN=alirc        ;ON=C2  ;VD=canopy albedo (near i.r.)                      ;VS=A*'//ncvp//'   ;VB=p0')
      PHYVAR3D1(alvsc,        'VN=alvsc        ;ON=C4  ;VD=canopy albedo (visible)                        ;VS=A*'//ncvp//'   ;VB=p0')
      PHYVAR3D1(bbi,          'VN=bbi          ;ON=C6  ;VD=Clapp and Hornberger hydraulic coefficient B   ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(cdh,          'VN=cdh          ;ON=G2  ;VD=bulk heat transfer coefficient                                    ;VB=v0')
      PHYVAR2D1(cdm,          'VN=cdm          ;ON=CM  ;VD=bulk momentum transfer coefficient                                ;VB=v0')
      PHYVAR3D1(clay,         'VN=clay         ;ON=CLAY;VD=percentage of clay in soil                     ;VS=A*'//ncg//'    ;VB=p1')
      PHYVAR2D1(cmai,         'VN=cmai         ;ON=CY  ;VD=instantaneous canopy mass                                         ;VB=p0')
      PHYVAR3D1(delzw,        'VN=delzw        ;ON=C7  ;VD=thickness of soil layers for water             ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(drain,        'VN=drain        ;ON=DR  ;VD=water drainage at bottom of soil layer                            ;VB=p0')
      PHYVAR2D1(drainaf,      'VN=drainaf      ;ON=O1  ;VD=accum. of base drainage                                           ;VB=p0')
      PHYVAR2D1(evapo,        'VN=evapo        ;ON=V5  ;VD=evapo./subl. rate from ground                                     ;VB=v0')
      PHYVAR3D1(fcanmx,       'VN=fcanmx       ;ON=Y2C ;VD=fract. coverage of vegetation classes          ;VS=A*'//ncvp//'   ;VB=p0')
      PHYVAR2D1(fcovc,        'VN=fcovc        ;ON=C3  ;VD=fractional coverage for canopy                                    ;VB=v0')
      PHYVAR2D1(fcovcs,       'VN=fcovcs       ;ON=S5  ;VD=fractional coverage for canopy+snow                               ;VB=v0')
      PHYVAR2D1(fcovg,        'VN=fcovg        ;ON=BG  ;VD=fractional coverage for bare ground                               ;VB=v0')
      PHYVAR2D1(fcovgs,       'VN=fcovgs       ;ON=S6  ;VD=fractional coverage for snow                                      ;VB=v0')
      PHYVAR2D1(firupaf,      'VN=firupaf      ;ON=N5  ;VD=acc. of soil surf. upward infrared flux                           ;VB=p0')
      PHYVAR2D1(flgg,         'VN=flgg         ;ON=F9  ;VD=LW radiation absorbed by ground                                   ;VB=v0')
      PHYVAR2D1(flgs,         'VN=flgs         ;ON=O3  ;VD=LW radiation absorbed by snow                                     ;VB=v0')
      PHYVAR2D1(flgv,         'VN=flgv         ;ON=F7  ;VD=LW radiation absorbed by canopy                                   ;VB=v0')
      PHYVAR2D1(fsgg,         'VN=fsgg         ;ON=F6  ;VD=SW radiation absorbed by ground                                   ;VB=v0')
      PHYVAR2D1(fsgs,         'VN=fsgs         ;ON=F5  ;VD=SW radiation absorbed by snow                                     ;VB=v0')
      PHYVAR2D1(fsgv,         'VN=fsgv         ;ON=F4  ;VD=SW radiation absorbed by canopy                                   ;VB=v0')
      PHYVAR2D1(fsnow,        'VN=fsnow        ;ON=5P  ;VD=diagnosed fractional snow coverage                                ;VB=v0')
      PHYVAR2D1(fsolupaf,     'VN=fsolupaf     ;ON=N6  ;VD=acc. of soil surf. upward solar flux                              ;VB=p0')
      PHYVAR2D1(grdhflx,      'VN=grdhflx      ;ON=GRHF;VD=Heat flux at soil surface [W/m^2]                                 ;VB=v0')
      PHYVAR2D1(grkfac,       'VN=grkfac       ;ON=GRKF;VD=WATROF par. for MESH code (grkfac)                                ;VB=v0')
      PHYVAR3D1(grksat,       'VN=grksat       ;ON=HT  ;VD=sat. soil hydraulic conductivity               ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(hcps,         'VN=hcps         ;ON=C9  ;VD=vol. heat capacity of soil                     ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(hevc,         'VN=hevc         ;ON=V2  ;VD=Latent heat flux from canopy                                      ;VB=v0')
      PHYVAR2D1(hevg,         'VN=hevg         ;ON=V4  ;VD=Latent heat flux from ground                                      ;VB=v0')
      PHYVAR2D1(hevs,         'VN=hevs         ;ON=V3  ;VD=Latent heat flux from snow                                        ;VB=v0')
      PHYVAR2D1(hfsc,         'VN=hfsc         ;ON=H4  ;VD=Sensible heat flux from canopy                                    ;VB=v0')
      PHYVAR2D1(hfsg,         'VN=hfsg         ;ON=H6  ;VD=Sensible heat flux from ground                                    ;VB=v0')
      PHYVAR2D1(hfss,         'VN=hfss         ;ON=H5  ;VD=Sensible heat flux from snow                                      ;VB=v0')
      PHYVAR2D1(hmfc,         'VN=hmfc         ;ON=G3  ;VD=Melting heat flux from canopy                                     ;VB=v0')
      PHYVAR3D1(hmfg,         'VN=hmfg         ;ON=O4  ;VD=Melting heat flux from ground                  ;VS=A*'//ncg//'    ;VB=v0')
      PHYVAR2D1(hmfn,         'VN=hmfn         ;ON=G4  ;VD=Melting heat flux from snow                                       ;VB=v0')
      PHYVAR3D1(htc,          'VN=htc          ;ON=M1C ;VD=Heat transfer through soil layers by perc.     ;VS=A*'//ncg//'    ;VB=v0')
      PHYVAR2D1(htcc,         'VN=htcc         ;ON=M2C ;VD=Heat transfer to veg. via pcp. interception                       ;VB=v0')
      PHYVAR2D1(htcs,         'VN=htcs         ;ON=Y3  ;VD=Heat transfer through snow via perc. and cond.                    ;VB=v0')
      PHYVAR2D1(huaircan,     'VN=huaircan     ;ON=SQAC;VD=specific humidity inside canopy                                   ;VB=p0')
      PHYVAR3D1(isoil,        'VN=isoil        ;ON=I2  ;VD=soil volumetric ice contents                   ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(iveg,         'VN=iveg         ;ON=SK  ;VD=snow stored on canopy                                             ;VB=p0')
      PHYVAR3D1(laimax,       'VN=laimax       ;ON=LX  ;VD=maximum leaf area index (LAI)                  ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(laimin,       'VN=laimin       ;ON=LN  ;VD=minimum leaf area index (LAI)                  ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(melts,        'VN=melts        ;ON=MLTS;VD=accum. snow melting (kg/m2)                                       ;VB=p0')
      PHYVAR2D1(mosfract,     'VN=mosfract     ;ON=MO  ;VD=amount of cover found in this mosaic tile                         ;VB=p0')
      PHYVAR3D1(orgm,         'VN=orgm         ;ON=Z7  ;VD=percentage of organic matter in soil           ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(overfl,       'VN=overfl       ;ON=RO  ;VD=overland runoff                                                   ;VB=v0')
      PHYVAR2D1(overflaf,     'VN=overflaf     ;ON=N0  ;VD=accum. of surface runoff                                          ;VB=p0')
      PHYVAR2D1(pcfc,         'VN=pcfc         ;ON=L5  ;VD=frozen precip. falling on canopy                                  ;VB=v0')
      PHYVAR2D1(pcpn,         'VN=pcpn         ;ON=PEG ;VD=precip incident on snow pack                                      ;VB=v0')
      PHYVAR2D1(pclc,         'VN=pclc         ;ON=P6  ;VD=liquid precip. falling on canopy                                  ;VB=v0')
      PHYVAR2D1(pcpg,         'VN=pcpg         ;ON=P7  ;VD=precip incident on ground                                         ;VB=v0')
      PHYVAR3D1(psiga,        'VN=psiga        ;ON=J3  ;VD=parameter psiga in stomatal resistance         ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(psigb,        'VN=psigb        ;ON=D4  ;VD=parameter psigb in stomatal resistance         ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(psisat,       'VN=psisat       ;ON=D5  ;VD=sat. soil water suction                        ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(psiwlt,       'VN=psiwlt       ;ON=D6  ;VD=soil water suction at wilting point            ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(qa50,         'VN=qa50         ;ON=D7  ;VD=parameter in stomatal conductance              ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(qfc,          'VN=qfc          ;ON=M5  ;VD=Water extract. from soil layers due to transp. ;VS=A*'//ncg//'    ;VB=v0')
      PHYVAR2D1(qfcf,         'VN=qfcf         ;ON=S1  ;VD=subl. rate of canopy frozen water                                 ;VB=v0')
      PHYVAR2D1(qfcl,         'VN=qfcl         ;ON=E2  ;VD=evapo. rate of canopy liq. water                                  ;VB=v0')
      PHYVAR2D1(qfg,          'VN=qfg          ;ON=E3  ;VD=evapo. rate from soil surface                                     ;VB=v0')
      PHYVAR2D1(qfn,          'VN=qfn          ;ON=S2  ;VD=subl. rate from snow cover                                        ;VB=v0')
      PHYVAR2D1(rib,          'VN=rib          ;ON=RIB ;VD=Bulk Richardson number [-10,5]                                    ;VB=v0')
      PHYVAR2D1(rofc,         'VN=rofc         ;ON=DC  ;VD=dripping from canopy                                              ;VB=v0')
      PHYVAR2D1(rofn,         'VN=rofn         ;ON=MS  ;VD=melting snow from snowpack                                        ;VB=v0')
      PHYVAR3D1(rootdp,       'VN=rootdp       ;ON=D2  ;VD=rooting soil depth                             ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(rovg,         'VN=rovg         ;ON=M7  ;VD=to be determined (rovg)                                           ;VB=v0')
      PHYVAR3D1(sand,         'VN=sand         ;ON=SAND;VD=percentage of sand in soil                     ;VS=A*'//ncg//'    ;VB=p1')
      PHYVAR2D1(snoden,       'VN=snoden       ;ON=DN  ;VD=snow density in kg/m3                                             ;VB=p1  ;MIN=0')
      PHYVAR2D1(soilcol,      'VN=soilcol      ;ON=SCOL;VD=soil color for albedo lookup table                                ;VB=p1')
      PHYVAR2D1(sdepth,       'VN=sdepth       ;ON=DPTH;VD=depth of soil water layer in CLASS                                ;VB=p1')
      PHYVAR2D1(snoal,        'VN=snoal        ;ON=I6  ;VD=albedo of snow                                                    ;VB=p0')
      PHYVAR2D1(snowsize,     'VN=snowsize     ;ON=SNGZ;VD=Snow grain size (for ISNOALB=1 option)  [m]                       ;VB=p0')
      PHYVAR3D1(stomr,        'VN=stomr        ;ON=RS  ;VD=minimum stomatal resistance                    ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(subflw,       'VN=subflw       ;ON=R6  ;VD=Interflow from sides of soil column                               ;VB=v0')
      PHYVAR2D1(taircan,      'VN=taircan      ;ON=STAC;VD=air temperature inside canopy                                     ;VB=p0')
      PHYVAR2D1(tbase,        'VN=tbase        ;ON=R2  ;VD=temp. at the base of soil water column                            ;VB=p0')
      PHYVAR2D1(tbasfl,       'VN=tbasfl       ;ON=STBF;VD=temperature of baseflow                                           ;VB=v0')
      PHYVAR3D1(tcs,          'VN=tcs          ;ON=H7  ;VD=thermal heat conductivity of soil              ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thfc,         'VN=thfc         ;ON=D9  ;VD=vol. water content at field capacity           ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thlmin,       'VN=thlmin       ;ON=E5  ;VD=minimun vol. water content                     ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thlrat,       'VN=thlrat       ;ON=E6  ;VD=soil retention capacity                        ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thlret,       'VN=thlret       ;ON=E7  ;VD=liq. water content behind wetting front        ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thpor,        'VN=thpor        ;ON=E8  ;VD=volumetric fraction of pores in soil           ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(tovrfl,       'VN=tovrfl       ;ON=STOF;VD=temperature of overland flow                                      ;VB=v0')
      PHYVAR2D1(tpond,        'VN=tpond        ;ON=Q4  ;VD=temperature of water lying on surface                             ;VB=p0')
      PHYVAR2D1(trunoff,      'VN=trunoff      ;ON=TROF;VD=temperature of runoff                                             ;VB=v0')
      PHYVAR2D1(tsno,         'VN=tsno         ;ON=TN  ;VD=ground snow temperature                                           ;VB=p0')
      PHYVAR3D1(tsoil,        'VN=tsoil        ;ON=I0  ;VD=surface and soil temperatures                  ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(tsubfl,       'VN=tsubfl       ;ON=STIF;VD=temperature of subsurface runoff                                  ;VB=v0')
      PHYVAR3D1(tsurfsa,      'VN=tsurfsa      ;ON=STSS;VD=surface tmp of class subareas                  ;VS=A*4            ;VB=p0')
      PHYVAR2D1(tveg,         'VN=tveg         ;ON=TE  ;VD=canopy temperature                                                ;VB=p0')
      PHYVAR2D1(veggro,       'VN=veggro       ;ON=GR  ;VD=canopy growth factor                                              ;VB=p0')
      PHYVAR3D1(vegma,        'VN=vegma        ;ON=MC  ;VD=standing mass of canopy                        ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(vpda,         'VN=vpda         ;ON=E9  ;VD=parameter vpda in stomatal resistance          ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(vpdb,         'VN=vpdb         ;ON=F1  ;VD=parameter vpdb in stomatal resistance          ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(wfcint,       'VN=wfcint       ;ON=WFCI;VD=WATROF par. for MESH code (wfcint)                                ;VB=v0')
      PHYVAR2D1(wflux,        'VN=wflux        ;ON=M8  ;VD=water flux from surface to atm.                                   ;VB=v0')
      PHYVAR2D1(wfluxaf,      'VN=wfluxaf      ;ON=N7  ;VD=acc. of soil surface upward water flux                            ;VB=p0')
      PHYVAR2D1(wfsurf,       'VN=wfsurf       ;ON=WFSF;VD=WATROF par. for MESH code (wfsurf)                                ;VB=v0')
      PHYVAR2D1(wsnow,        'VN=wsnow        ;ON=I4  ;VD=water in the snow pack                                            ;VB=p0')
      PHYVAR3D1(wsoil,        'VN=wsoil        ;ON=I1  ;VD=soil volumetric water contents                 ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(wtrc,         'VN=wtrc         ;ON=W6C ;VD=water lost from canopy(senescence)                                ;VB=v0')
      PHYVAR2D1(wtrg,         'VN=wtrg         ;ON=W8C ;VD=water lost to surface                                             ;VB=v0')
      PHYVAR2D1(wtrs,         'VN=wtrs         ;ON=W7C ;VD=water lost to snowpack                                            ;VB=v0')
      PHYVAR2D1(wveg,         'VN=wveg         ;ON=I3  ;VD=water retained on the vegetation                                  ;VB=p0')
      PHYVAR2D1(xdrain,       'VN=xdrain       ;ON=L9  ;VD=drainage factor in CLASS                                          ;VB=p0')
      PHYVAR2D1(xslope,       'VN=xslope       ;ON=SFIS;VD=mosaic slope                                                      ;VB=p0')
      PHYVAR3D1(zbotw,        'VN=zbotw        ;ON=G0  ;VD=soil layer depths for water                    ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(zoln,         'VN=zoln         ;ON=X9  ;VD=roughness length for each veg. class           ;VS=A*'//ncvp//'   ;VB=p0')
      PHYVAR2D1(zpond,        'VN=zpond        ;ON=M9  ;VD=height of water lying on surface                                  ;VB=p0')

      ! New variables for running the alternate ground water scheme
      PHYVAR2D1(anis,         'VN=anis         ;ON=ANIS;VD=anisotropy                                                        ;VB=p0')
      PHYVAR2D1(are,          'VN=are          ;ON=ARE ;VD=grid point area (km2)                                             ;VB=p0')
      PHYVAR2D1(excw,         'VN=excw         ;ON=EXCW;VD=water depth exceeding ground water capacity                       ;VB=p0')
      PHYVAR2D1(lbedr,        'VN=lbedr        ;ON=IBED;VD=soil layer containing water table                                 ;VB=p0')
      PHYVAR2D1(leggw,        'VN=leggw        ;ON=LEG ;VD=canal length                                                      ;VB=p0')
      PHYVAR2D1(slpgw,        'VN=slpgw        ;ON=SLP ;VD=land surface slope                                                ;VB=p0')
      PHYVAR2D1(totw,         'VN=totw         ;ON=WATT;VD=ground water storage                                              ;VB=p0')
      PHYVAR2D1(wtnew,        'VN=wtnew        ;ON=WT2 ;VD=updated ground water storage                                      ;VB=p0')

! Parameters for CLASS 3.6
      PHYVAR2D1(algwet,       'VN=algwet       ;ON=ALWT;VD=albedo of wet soil                                                ;VB=p0')
      PHYVAR2D1(algdry,       'VN=algdry       ;ON=ALDY;VD=albedo of dry soil                                                ;VB=p0')

      ! Parameters for CTEM
      PHYVAR3D1(ailc,         'VN=ailc         ;ON=AILC;VD=lumped LAI for 4 CLASS PFTs                    ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(ailcb,        'VN=ailcb        ;ON=AILB;VD=brown LAI                                      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(ailcg,        'VN=ailcg        ;ON=AILG;VD=green LAI                                      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(alirctm,      'VN=alirctm      ;ON=ALIC;VD=canopy albedo (near i.r.) from CTEM            ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(allwacc,      'VN=allwacc      ;ON=LAAC;VD=longwave albedo acc. for cur. day                                 ;VB=p0')
      PHYVAR2D1(alswacc,      'VN=alswacc      ;ON=SAAC;VD=shortwave albedo acc. for cur. day                                ;VB=p0')
      PHYVAR3D1(alvsctm,      'VN=alvsctm      ;ON=ALIS;VD=canopy albedo (visible) from CTEM              ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(ancgvgac,     'VN=ancgvgac     ;ON=PHGA;VD=daily accum. of ancgveg                        ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(ancsvgac,     'VN=ancsvgac     ;ON=PHSA;VD=daily accum. of ancsveg                        ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(anndefct,     'VN=anndefct     ;ON=AWD ;VD=an. water deficit (mm)                                            ;VB=p0')
      PHYVAR2D1(annpcp,       'VN=annpcp       ;ON=APR ;VD=an. prec. (mm)                                                    ;VB=p0')
      PHYVAR2D1(annsrpls,     'VN=annsrpls     ;ON=AWS ;VD=an. water surplus (mm)                                            ;VB=p0')
      PHYVAR2D1(anpcpcur,     'VN=anpcpcur     ;ON=APRY;VD=annual prec. for cur. year (mm)                                   ;VB=p0')
      PHYVAR2D1(anpecur,      'VN=anpecur      ;ON=AEVY;VD=an. pot. evap. for cur. year (mm)                                 ;VB=p0')
      PHYVAR2D1(anpotevp,     'VN=anpotevp     ;ON=AEV ;VD=an. pot. evap. (mm)                                               ;VB=p0')
      PHYVAR2D1(aridity,      'VN=aridity      ;ON=ARI ;VD=aridity index                                                     ;VB=p0')
      PHYVAR3D1(bleafmas,     'VN=bleafmas     ;ON=CBLF;VD=brown leaf mass                                ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(bmasveg,      'VN=bmasveg      ;ON=BMV ;VD=total (gleaf + stem + root) biomass            ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(burnvegf,     'VN=burnvegf     ;ON=BVF ;VD=burned area fraction                           ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(cfluxcg,      'VN=cfluxcg      ;ON=CFCG;VD=aerodynamic conductance over ground                               ;VB=p0')
      PHYVAR2D1(cfluxcs,      'VN=cfluxcs      ;ON=CFCS;VD=aerodynamic conductance over snow                                 ;VB=p0')
      PHYVAR3D1(cmasvegc,     'VN=cmasvegc     ;ON=CMVC;VD=total canopy mass from CTEM                    ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR2D1(co2conc,      'VN=co2conc      ;ON=CO2C;VD=atmos. CO2 conc.                                                  ;VB=p0')
      PHYVAR3D1(co2i1cg,      'VN=co2i1cg      ;ON=C1G ;VD=intercellular CO2 conc. (ground, sunlit)       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(co2i1cs,      'VN=co2i1cs      ;ON=C1S ;VD=intercellular CO2 conc. (snow, sunlit)         ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(co2i2cg,      'VN=co2i2cg      ;ON=C2G ;VD=intercellular CO2 conc. (ground, shaded)       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(co2i2cs,      'VN=co2i2cs      ;ON=C2S ;VD=intercellular CO2 conc. (snow, shaded)         ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(colddayr,     'VN=colddayr     ;ON=CLDD;VD=cold days counter                              ;VS=A*2            ;VB=p0')
      PHYVAR2D1(defctcur,     'VN=defctcur     ;ON=WDM ;VD=water deficit for cur. month                                      ;VB=p0')
      PHYVAR2D1(defctmon,     'VN=defctmon     ;ON=MWD ;VD=nb months with water deficit                                      ;VB=p0')
      PHYVAR2D1(defmnr,       'VN=defmnr       ;ON=MWDY;VD=nb months with water deficit for cur. yr                          ;VB=p0')
      PHYVAR2D1(dftcuryr,     'VN=dftcuryr     ;ON=AWDY;VD=water deficit for cur. year                                       ;VB=p0')
      PHYVAR2D1(dryslen,      'VN=dryslen      ;ON=DSL ;VD=dry season length (months)                                        ;VB=p0')
      PHYVAR3D1(dvdfcan,      'VN=dvdfcan      ;ON=DVDF;VD=CTEM PFT fractions within CLASS PFTs           ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(extnprob,     'VN=extnprob     ;ON=XTP ;VD=fire extinguishing probability                                    ;VB=p0')
      PHYVAR3D1(fcancmx,      'VN=fcancmx      ;ON=FCAN; '//'VD=max. fract. coverage of each CTEM PFT     ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(flhrloss,     'VN=flhrloss     ;ON=FHL ;VD=fall or harvest loss                           ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(flinacc,      'VN=flinacc      ;ON=FIAC;VD=IR downward energy flux acc. for cur. day                         ;VB=p0')
      PHYVAR2D1(flutacc,      'VN=flutacc      ;ON=FTAC;VD=IR upward energy flux acc. for cur. day                           ;VB=p0')
      PHYVAR2D1(fsinacc,      'VN=fsinacc      ;ON=FBAC;VD=VIS downward flux acc. for the day                                ;VB=p0')
      PHYVAR2D1(fsnowacc,     'VN=fsnowacc     ;ON=FSNA;VD=daily accum. snow fraction                                        ;VB=p0')
      PHYVAR2D1(gavglai,      'VN=gavglai      ;ON=GLAI;VD=grid averaged green LAI                                           ;VB=p0')
      PHYVAR2D1(gavgltms,     'VN=gavgltms     ;ON=GLF ;VD=grid averaged litter mass                                         ;VB=p0')
      PHYVAR2D1(gavgscms,     'VN=gavgscms     ;ON=GSTM;VD=grid averaged soil carbon mass                                    ;VB=p0')
      PHYVAR2D1(gdd5,         'VN=gdd5         ;ON=GD5 ;VD=growing degree days over 5c                                       ;VB=p0')
      PHYVAR2D1(gdd5cur,      'VN=gdd5cur      ;ON=GD5Y;VD=growing degree days over 5c for cur. yr                           ;VB=p0')
      PHYVAR3D1(geremort,     'VN=geremort     ;ON=GMRT;VD=growth related mortality (1/day)               ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(gleafmas,     'VN=gleafmas     ;ON=CGLF;VD=green leaf mass                                ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(grwtheff,     'VN=grwtheff     ;ON=GREF;VD=growth efficiency                              ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(intrmort,     'VN=intrmort     ;ON=IMRT;VD=intrinsic (age related) mortality (1/day)      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(lambda,       'VN=lambda       ;ON=LMBD;VD=npp fraction used for spatial expansion        ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(lfstatur,     'VN=lfstatur     ;ON=HLST;VD=leaf phenology status                          ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(litrmass,     'VN=litrmass     ;ON=CLTR;VD=litter mass                                    ;VS=A*'//niccp//'  ;VB=p0')
      PHYVAR3D1(lyrotmas,     'VN=lyrotmas     ;ON=LYRM;VD=root mass at the end of last year              ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(lystmmas,     'VN=lystmmas     ;ON=LYSM;VD=stem mass at the end of last year              ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(mlightng,     'VN=mlightng     ;ON=MLNT;VD=mean monthly lightning frequency               ;VS=A*12           ;VB=p0')
      PHYVAR3D1(nfcancmx,     'VN=nfcancmx     ;ON=NY2C;VD=fcancmx for next year                          ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(nppveg,       'VN=nppveg       ;ON=NPPV;VD=net primary productivity                       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(paic,         'VN=paic         ;ON=PAIC;VD=PAI for CLASS PFTs                             ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(pandayr,      'VN=pandayr      ;ON=HPAD;VD=days with positive net photosynthesis          ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(pfcancmx,     'VN=pfcancmx     ;ON=PY2C;VD=fcancmx from previous year                     ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(pftexistr,    'VN=pftexistr    ;ON=PFTE;VD=array indic. pfts exist (1) or not (-1)        ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(pgleafmass,   'VN=pgleafmass   ;ON=PGLM;VD=root mass (?) from prev. step                  ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(prbfrhuc,     'VN=prbfrhuc     ;ON=PFHC;VD=prob. of fire due to human causes                                 ;VB=p0')
      PHYVAR2D1(preacc,       'VN=preacc       ;ON=PRAC;VD=accum. of precip. for the day                                     ;VB=p0')
      PHYVAR3D1(pstemmass,    'VN=pstemmass    ;ON=PSTM;VD=stem mass from prev. step                      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(rmatc,        'VN=rmatc        ;ON=RMA ;VD=layer root fraction (CLASS PFTs)               ;VS=A*'//ncvxcg//' ;VB=p0')
      PHYVAR3D1(rmatctem,     'VN=rmatctem     ;ON=RMAC;VD=layer root fraction (CTEM PFTs)                ;VS=A*'//niccxcg//';VB=p0')
      PHYVAR3D1(rmlcgvga,     'VN=rmlcgvga     ;ON=LRGA;VD=daily accum. of rmlcgveg                       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(rmlcsvga,     'VN=rmlcsvga     ;ON=LRSA;VD=daily accum. of rmlcsveg                       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(rootdpth,     'VN=rootdpth     ;ON=RTDP;VD=rooting depth                                  ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(rootmass,     'VN=rootmass     ;ON=CROT;VD=root mass                                      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(rothrlos,     'VN=rothrlos     ;ON=RHL ;VD=root loss due to harvest                       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(slai,         'VN=slai         ;ON=SLAI;VD=imaginary LAI for phenology purposes           ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(slaic,        'VN=slaic        ;ON=SLC ;VD=storage LAI for use within CLASS               ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(soilcmas,     'VN=soilcmas     ;ON=CSOC;VD=soil carbon mass                               ;VS=A*'//niccp//'  ;VB=p0')
      PHYVAR2D1(srpcuryr,     'VN=srpcuryr     ;ON=AWSY;VD=water surplus for cur. yr                                         ;VB=p0')
      PHYVAR2D1(srplscur,     'VN=srplscur     ;ON=WSM ;VD=water surplus for cur. month                                      ;VB=p0')
      PHYVAR2D1(srplsmon,     'VN=srplsmon     ;ON=MWS ;VD=nb months with water surplus                                      ;VB=p0')
      PHYVAR3D1(stemmass,     'VN=stemmass     ;ON=CSTM;VD=stem mass                                      ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR3D1(stmhrlos,     'VN=stmhrlos     ;ON=SHL ;VD=stem loss due to harvest                       ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(surmnr,       'VN=surmnr       ;ON=MWSY;VD=nb months with water surplus for cur. yr                          ;VB=p0')
      PHYVAR2D1(taaccgat,     'VN=taaccgat     ;ON=TAAC;VD=daily accum. of air temperature                                   ;VB=p0')
      PHYVAR3D1(tbaraccgat,   'VN=tbaraccgat   ;ON=TBAC;VD=daily accum. of soil temperature               ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(tbarcacc,     'VN=tbarcacc     ;ON=TCA ;VD=soil temp. acc. (canopy over ground)           ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(tbarcsacc,    'VN=tbarcsacc    ;ON=TCSA;VD=soil temp. acc. (canopy over snow)             ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(tbargacc,     'VN=tbargacc     ;ON=TGA ;VD=soil temp. acc. (bare ground)                  ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(tbargsacc,    'VN=tbargsacc    ;ON=TGSA;VD=soil temp. acc. (snow over ground)             ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR2D1(tcanoaccgat,  'VN=tcanoaccgat  ;ON=TCOA;VD=canopy temp. acc. over ground                                     ;VB=p0')
      PHYVAR2D1(tcansacc,     'VN=tcansacc     ;ON=TSA ;VD=canopy temp. acc. over snow                                       ;VB=p0')
      PHYVAR2D1(tcoldm,       'VN=tcoldm       ;ON=CMT ;VD=temperature of coldest month (c)                                  ;VB=p0')
      PHYVAR2D1(tcurm,        'VN=tcurm        ;ON=TCM ;VD=temperature of current month (c)                                  ;VB=p0')
      PHYVAR3D1(thicecacc,    'VN=thicecacc    ;ON=THIA;VD=daily acc. canopy frozen water                 ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thliqcacc,    'VN=thliqcacc    ;ON=TLCA;VD=daily acc. canopy liquid water                 ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(thliqgacc,    'VN=thliqgacc    ;ON=TLGA;VD=daily acc. liquid water on ground              ;VS=A*'//ncg//'    ;VB=p0')
      PHYVAR3D1(tmonthb,      'VN=tmonthb      ;ON=T12 ;VD=monthly temperatures (c)                       ;VS=A*12           ;VB=p0')
      PHYVAR3D1(todfrac,      'VN=todfrac      ;ON=TODF;VD=max coverage of CTEM PFTs at end of day        ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(twarmm,       'VN=twarmm       ;ON=WMT ;VD=temperature of warmest month (c)                                  ;VB=p0')
      PHYVAR3D1(tymaxlai,     'VN=tymaxlai     ;ON=TYML;VD=max LAI for this year)                         ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(uvaccgat,     'VN=uvaccgat     ;ON=UVAC;VD=daily acc. U wind speed                                           ;VB=p0')
      PHYVAR3D1(veghght,      'VN=veghght      ;ON=VGHG;VD=vegetation height                              ;VS=A*'//nicc//'   ;VB=p0')
      PHYVAR2D1(vgbiomas,     'VN=vgbiomas     ;ON=VGBM;VD=grid averaged vegetation biomass                                  ;VB=p0')
      PHYVAR2D1(vvaccgat,     'VN=vvaccgat     ;ON=VVAC;VD=daily acc. V wind speed                                           ;VB=p0')
      PHYVAR3D1(wdmindex,     'VN=wdmindex     ;ON=WDI ;VD=array indic. wet (1) or dry (-1) month         ;VS=A*12           ;VB=p0')
      PHYVAR3D1(zolnc,        'VN=zolnc        ;ON=Z0C ;VD=roughness length for CLASS PFTs                ;VS=A*'//ncv//'    ;VB=p0')
      PHYVAR3D1(afrleaf,      'VN=afrleaf      ;ON=AFRL;VD=leaf allocation fraction                       ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(afrroot,      'VN=afrroot      ;ON=AFRR;VD=root allocation fraction                       ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(afrstem,      'VN=afrstem      ;ON=AFRS;VD=stem allocation fraction                       ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(autores,      'VN=autores      ;ON=AUTR;VD=grid avg. autotrophic resp.                                       ;VB=v0')
      PHYVAR3D1(autresveg,    'VN=autresveg    ;ON=ATRV;VD=autotrophic respiration                        ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(burnarea,     'VN=burnarea     ;ON=BURN;VD=burned area                                                       ;VB=v0')
      PHYVAR3D1(colrate,      'VN=colrate      ;ON=COL ;VD=colonisation rate (1/day)                      ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(dstcemls,     'VN=dstcemls     ;ON=DST ;VD=co2 emission losses from veg. disturbance                         ;VB=v0')
      PHYVAR2D1(dstcemls3,    'VN=dstcemls3    ;ON=DST3;VD=co2 emission losses from litter pool dist.                        ;VB=v0')
      PHYVAR2D1(gpp,          'VN=gpp          ;ON=GPP ;VD=grid avg. gross primary productivity                              ;VB=v0')
      PHYVAR3D1(gppveg,       'VN=gppveg       ;ON=GPPV;VD=gross primary productivity                     ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(grclarea,     'VN=grclarea     ;ON=GRCL;VD=grid cell area                                                    ;VB=v0')
      PHYVAR3D1(hetresveg,    'VN=hetresveg    ;ON=HTRV;VD=heterotrophic respiration                      ;VS=A*'//niccp//'  ;VB=v0')
      PHYVAR2D1(hetrores,     'VN=hetrores     ;ON=HTRR;VD=grid avg. heterotrophic resp.                                     ;VB=v0')
      PHYVAR2D1(humiftrs,     'VN=humiftrs     ;ON=HUMF;VD=humidified litter transfer to soil C pool                         ;VB=v0')
      PHYVAR3D1(leaflitr,     'VN=leaflitr     ;ON=LLIT;VD=leaf litter fall not caused by fire or mort.   ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(litrfall,     'VN=litrfall     ;ON=LITF;VD=total litter fall                                                 ;VB=v0')
      PHYVAR2D1(litres,       'VN=litres       ;ON=LITR;VD=grid avg. litter respiration                                      ;VB=v0')
      PHYVAR3D1(litresveg,    'VN=litresveg    ;ON=LTRV;VD=litter respiration                             ;VS=A*'//niccp//'  ;VB=v0')
      PHYVAR2D1(lucemcom,     'VN=lucemcom     ;ON=LUCE;VD=land use change combustion emission losses                        ;VB=v0')
      PHYVAR2D1(lucltrin,     'VN=lucltrin     ;ON=LUCL;VD=land use change litter pool inputs                                ;VB=v0')
      PHYVAR2D1(lucsocin,     'VN=lucsocin     ;ON=LUCS;VD=land use change soil C pool inputs                                ;VB=v0')
      PHYVAR3D1(ltstatus,     'VN=ltstatus     ;ON=LTST;VD=light status                                   ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(mortrate,     'VN=mortrate     ;ON=MORT;VD=mortality rate (1/day)                         ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(nbp,          'VN=nbp          ;ON=NBP ;VD=grid averaged net biome productivity                              ;VB=v0')
      PHYVAR3D1(nbpveg,       'VN=nbpveg       ;ON=NBPV;VD=net biome productivity                         ;VS=A*'//niccp//'  ;VB=v0')
      PHYVAR2D1(nep,          'VN=nep          ;ON=NEP ;VD=grid averaged net ecosystem productivity                          ;VB=v0')
      PHYVAR3D1(nepveg,       'VN=nepveg       ;ON=NEPV;VD=net ecosystem productivity                     ;VS=A*'//niccp//'  ;VB=v0')
      PHYVAR2D1(npp,          'VN=npp          ;ON=NPP ;VD=grid averaved net primary productivity                            ;VB=v0')
      PHYVAR2D1(probfire,     'VN=probfire     ;ON=PRBF;VD=probability of fire                                               ;VB=v0')
      PHYVAR2D1(rg,           'VN=rg           ;ON=RGC ;VD=grid averaged growth respiration                                  ;VB=v0')
      PHYVAR3D1(rgveg,        'VN=rgveg        ;ON=RGCV;VD=growth respiration                             ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(rm,           'VN=rm           ;ON=RMC ;VD=maintenance respiration                                           ;VB=v0')
      PHYVAR2D1(rml,          'VN=rml          ;ON=RML ;VD=grid avg. leaf maintenance resp.                                  ;VB=v0')
      PHYVAR3D1(rmlvegacc,    'VN=rmlvegacc    ;ON=RMLA;VD=leaf maintenance respiration                   ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(rmr,          'VN=rmr          ;ON=RMR ;VD=grid avg. root maintenance resp.                                  ;VB=v0')
      PHYVAR3D1(rmrveg,       'VN=rmrveg       ;ON=RMRV;VD=root maintenance resp.                         ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(rms,          'VN=rms          ;ON=RMS ;VD=grid avg. stem maintenance resp.                                  ;VB=v0')
      PHYVAR3D1(rmsveg,       'VN=rmsveg       ;ON=RMSV;VD=stem maintenance resp.                         ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(roottemp,     'VN=roottemp     ;ON=RTC ;VD=root temperature                               ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR2D1(socres,       'VN=socres       ;ON=SOCR;VD=grid avg. soil carbon resp.                                       ;VB=v0')
      PHYVAR3D1(socresveg,    'VN=socresveg    ;ON=SCRV;VD=soil carbon respiration                        ;VS=A*'//niccp//'  ;VB=v0')
      PHYVAR2D1(soilresp,     'VN=soilresp     ;ON=SOLR;VD=soil respiration                                                  ;VB=v0')
      PHYVAR3D1(tltrleaf,     'VN=tltrleaf     ;ON=TLTL;VD=leaf litter fall                               ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(tltrroot,     'VN=tltrroot     ;ON=TLTR;VD=root litter fall                               ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(tltrstem,     'VN=tltrstem     ;ON=TLTS;VD=stem litter fall                               ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(vgbiomas_veg, 'VN=vgbiomas_veg ;ON=VGBG;VD=vegetation biomass                             ;VS=A*'//nicc//'   ;VB=v0')
      PHYVAR3D1(wtstatus,     'VN=wtstatus     ;ON=WTST;VD=soil water status                              ;VS=A*'//nicc//'   ;VB=v0')
   endif IF_CLASS


   IF_LAKES: if (schmlake /= 'NIL') then
      PHYVAR2D1(frv_li,       'VN=frv_li       ;ON=LIFV;VD=Lake ice surface friction velocity                                ;VB=p0')
      PHYVAR2D1(frv_lw,       'VN=frv_lw       ;ON=LWFV;VD=Lake water surface friction velocity                              ;VB=p0')
      PHYVAR2D1(lakedepth,    'VN=lakedepth    ;ON=LDEP;VD=Lake depth                                                        ;VB=p1')
      PHYVAR2D1(lakefice,     'VN=lakefice     ;ON=LIFR;VD=Lake ice fraction                                                 ;VB=p0')
      PHYVAR2D1(lakehice,     'VN=lakehice     ;ON=LITH;VD=Lake ice thickness                                                ;VB=p0')
      PHYVAR2D1(laketice,     'VN=laketice     ;ON=LITP;VD=Lake ice surface temperature                                      ;VB=p0')
      PHYVAR2D1(laketransp,   'VN=laketransp   ;ON=LTRN;VD=Lake water transparency                                           ;VB=p0')
      PHYVAR2D1(lakehml,      'VN=lakehml      ;ON=LMLD;VD=Lake mixed layer thickness                                        ;VB=p0')
      if (schmlake == 'FLAKE') then
         PHYVAR2D1(lakect,    'VN=lakect       ;ON=LSF ;VD=Lake shape factor (thermocline) in FLake                          ;VB=p0')
         PHYVAR2D1(laketbot,  'VN=laketbot     ;ON=LBTP;VD=Lake bottom temperature in FLake                                  ;VB=p0')
         PHYVAR2D1(laketmnw,  'VN=laketmnw     ;ON=LWTP;VD=Lake water average temperature in FLake                           ;VB=p0')
         PHYVAR2D1(laketwml,  'VN=laketwml     ;ON=LMLT;VD=Lake mixed layer temperature in FLake                             ;VB=p0')
      endif
   endif IF_LAKES

   if (schmurb /= 'TEB') return

   PHYVAR2D1(alb_road,     'VN=alb_road     ;ON=ALRD;VD=road albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_roaden,   'VN=alb_roaden   ;ON=TB9 ;VD=road albedo (E)                                           ;VB=e1; IN=ALRD;MIN=0')
   PHYVAR2D1(alb_roof,     'VN=alb_roof     ;ON=ALRF;VD=roof albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_roofen,   'VN=alb_roofen   ;ON=TB10;VD=roof albedo (E)                                           ;VB=e1; IN=ALRF;MIN=0')
   PHYVAR2D1(alscatw,      'VN=alscatw      ;ON=ALSC;VD=Town albedo for scattered solar radiation                 ;VB=p0')
   PHYVAR2D1(alb_wall,     'VN=alb_wall     ;ON=ALWL;VD=wall albedo                                               ;VB=p0         ;MIN=0')
   PHYVAR2D1(alb_wallen,   'VN=alb_wallen   ;ON=TB11;VD=wall albedo (E)                                           ;VB=e1; IN=ALWL;MIN=0')
   PHYVAR2D1(azim,         'VN=azim         ;ON=AZIM;VD=solar azimuthal angle                                     ;VB=p0')
   PHYVAR2D1(bld,          'VN=bld          ;ON=BLDF;VD=building fraction                                         ;VB=p0         ;MIN=0')
   PHYVAR2D1(blden,        'VN=blden        ;ON=TB2 ;VD=building fraction (E)                                     ;VB=e1; IN=BLDF;MIN=0')
   PHYVAR2D1(bld_height,   'VN=bld_height   ;ON=BLDH;VD=building height                                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(bld_heighten, 'VN=bld_heighten ;ON=TB3 ;VD=building height (E)                                       ;VB=e1; IN=BLDH;MIN=0')
   PHYVAR2D1(can_hw_ratio, 'VN=can_hw_ratio ;ON=ASPC;VD=aspect ratio of the street                                ;VB=p0         ;MIN=0')
   PHYVAR3D1(d_road,       'VN=d_road       ;ON=DPRD;VD=depth of the road layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_roaden,     'VN=d_roaden     ;ON=TB12;VD=depth of the road layers (E)                   ;VS=A*3  ;VB=e1; IN=DPRD;')
   PHYVAR3D1(d_roof,       'VN=d_roof       ;ON=DPRF;VD=depth of the roof layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_roofen,     'VN=d_roofen     ;ON=TB13;VD=depth of the roof layers (E)                   ;VS=A*3  ;VB=e1; IN=DPRF;')
   PHYVAR3D1(d_wall,       'VN=d_wall       ;ON=DPWL;VD=depth of the wall layers                       ;VS=A*3  ;VB=p0')
   PHYVAR3D1(d_wallen,     'VN=d_wallen     ;ON=TB14;VD=depth of the wall layers (E)                   ;VS=A*3  ;VB=e1; IN=DPWL;')
   PHYVAR2D1(emis_road,    'VN=emis_road    ;ON=EMRD;VD=road emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_roaden,  'VN=emis_roaden  ;ON=TB15;VD=road emissivity (E)                                       ;VB=e1; IN=EMRD;')
   PHYVAR2D1(emis_roof,    'VN=emis_roof    ;ON=EMRF;VD=roof emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_roofen,  'VN=emis_roofen  ;ON=TB16;VD=roof emissivity (E)                                       ;VB=e1; IN=EMRF;')
   PHYVAR2D1(emis_wall,    'VN=emis_wall    ;ON=EMWL;VD=wall emissivity                                           ;VB=p0')
   PHYVAR2D1(emis_wallen,  'VN=emis_wallen  ;ON=TB17;VD=wall emissivity (E)                                       ;VB=e1; IN=EMWL;')
   PHYVAR2D1(emtw,         'VN=emtw         ;ON=EMTW;VD=Town emissivity                                           ;VB=p0')
   PHYVAR2D1(g_road,       'VN=g_road       ;ON=QGRD;VD=storage heat flux for road                                ;VB=p0')
   PHYVAR2D1(g_roof,       'VN=g_roof       ;ON=QGRF;VD=storage heat flux for roof                                ;VB=p0')
   PHYVAR2D1(g_town,       'VN=g_town       ;ON=QGTW;VD=storage heat flux for town                                ;VB=p0')
   PHYVAR2D1(g_wall,       'VN=g_wall       ;ON=QGWL;VD=storage heat flux for wall                                ;VB=p0')
   PHYVAR2D1(h_industry,   'VN=h_industry   ;ON=QHIN;VD=sensible heat flux from industry                          ;VB=p0')
   PHYVAR2D1(h_industryen, 'VN=h_industryen ;ON=TB18;VD=sensible heat flux from industry (E)                      ;VB=e1; IN=QHIN;')
   PHYVAR2D1(h_road,       'VN=h_road       ;ON=QHRD;VD=sensible heat flux over road                              ;VB=p0')
   PHYVAR2D1(h_roof,       'VN=h_roof       ;ON=QHRF;VD=sensible heat flux over roof                              ;VB=p0')
   PHYVAR2D1(h_town,       'VN=h_town       ;ON=QHTW;VD=sensible heat flux over town                              ;VB=p0')
   PHYVAR2D1(h_traffic,    'VN=h_traffic    ;ON=QHTR;VD=sensible heat flux from traffic                           ;VB=p0')
   PHYVAR2D1(h_trafficen,  'VN=h_trafficen  ;ON=TB28;VD=sensible heat flux from traffic (E)                       ;VB=e1; IN=QHTR;')
   PHYVAR2D1(h_wall,       'VN=h_wall       ;ON=QHWL;VD=sensible heat flux over wall                              ;VB=p0')
   PHYVAR3D1(hc_road,      'VN=hc_road      ;ON=HCRD;VD=road heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_roaden,    'VN=hc_roaden    ;ON=TB19;VD=road heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCRD;')
   PHYVAR3D1(hc_roof,      'VN=hc_roof      ;ON=HCRF;VD=roof heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_roofen,    'VN=hc_roofen    ;ON=TB20;VD=roof heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCRF;')
   PHYVAR3D1(hc_wall,      'VN=hc_wall      ;ON=HCWL;VD=wall heat capacities                           ;VS=A*3  ;VB=p0')
   PHYVAR3D1(hc_wallen,    'VN=hc_wallen    ;ON=TB21;VD=wall heat capacities (E)                       ;VS=A*3  ;VB=e1; IN=HCWL;')
   PHYVAR2D1(le_industry,  'VN=le_industry  ;ON=QEIN;VD=latent heat flux from industry                            ;VB=p0')
   PHYVAR2D1(le_industryen,'VN=le_industryen;ON=TB27;VD=latent heat flux from industry (E)                        ;VB=e1; IN=QEIN;')
   PHYVAR2D1(le_road,      'VN=le_road      ;ON=QERD;VD=latent heat flux over road                                ;VB=p0')
   PHYVAR2D1(le_roof,      'VN=le_roof      ;ON=QERF;VD=latent heat flux over roof                                ;VB=p0')
   PHYVAR2D1(le_town,      'VN=le_town      ;ON=QETW;VD=latent heat flux over town                                ;VB=p0')
   PHYVAR2D1(le_traffic,   'VN=le_traffic   ;ON=QETR;VD=latent heat flux from traffic                             ;VB=p0')
   PHYVAR2D1(le_trafficen, 'VN=le_trafficen ;ON=TB26;VD=latent heat flux from traffic (E)                         ;VB=e1; IN=QETR;')
   PHYVAR2D1(le_wall,      'VN=le_wall      ;ON=QEWL;VD=latent heat flux over wall                                ;VB=p0')
   PHYVAR2D1(nat,          'VN=nat          ;ON=NATF;VD=natural surface fraction in urban area                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(naten,        'VN=naten        ;ON=TB1 ;VD=natural surface fraction in urban area (E)                ;VB=e1; IN=NATF;MIN=0')
   PHYVAR2D1(pav,          'VN=pav          ;ON=PAVF;VD=impervious fraction (road)                                ;VB=p0         ;MIN=0')
   PHYVAR2D1(paven,        'VN=paven        ;ON=TB4 ;VD=impervious fraction (road) (E)                            ;VB=e1; IN=PAVF;MIN=0')
   PHYVAR2D1(q_canyon,     'VN=q_canyon     ;ON=QCAN;VD=specific humi inside the canyon                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(q_canyonen,   'VN=q_canyonen   ;ON=2XEN; VD=canyon air humi (E)                                      ;VB=e1 ;IN=QCAN;MIN=0')
   PHYVAR2D1(rn_road,      'VN=rn_road      ;ON=RNRD;VD=Net radiation over road                                   ;VB=p0')
   PHYVAR2D1(rn_roof,      'VN=rn_roof      ;ON=RNRF;VD=Net radiation over roof                                   ;VB=p0')
   PHYVAR2D1(rn_town,      'VN=rn_town      ;ON=RNTW;VD=Net radiation over town                                   ;VB=p0')
   PHYVAR2D1(rn_wall,      'VN=rn_wall      ;ON=RNWL;VD=Net radiation over wall                                   ;VB=p0')
   PHYVAR2D1(sroad_alb,    'VN=sroad_alb    ;ON=SARD;VD=snow albedo for roads                                     ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_alben,  'VN=sroad_alben  ;ON=9ZEN; VD=snow albedo for roads (E)                                ;VB=e1 ;IN=SARD;MIN=0')
   PHYVAR2D1(sroad_emis,   'VN=sroad_emis   ;ON=SERD;VD=snow emissivity for roads                                 ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_emisen, 'VN=sroad_emisen ;ON=1MEN; VD=snow emmissivity for roads (E)                           ;VB=e1 ;IN=SERD;MIN=0')
   PHYVAR2D1(sroad_rho,    'VN=sroad_rho    ;ON=SDRD;VD=snow density for roads                                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroad_rhoen,  'VN=sroad_rhoen  ;ON=8ZEN; VD=snow density for roads (E)                               ;VB=e1 ;IN=SDRD;MIN=0')
   PHYVAR2D1(sroad_scheme, 'VN=sroad_scheme ;ON=SCRD;VD=snow scheme for roads                                     ;VB=p0')
   PHYVAR2D1(sroad_t,      'VN=sroad_t      ;ON=STRD;VD=snow temperature for roads                                ;VB=p0')
   PHYVAR2D1(sroad_ten,    'VN=sroad_ten    ;ON=7ZEN; VD=snow temp for roads (E)                                  ;VB=e1 ;IN=STRD;')
   PHYVAR2D1(sroad_ts,     'VN=sroad_ts     ;ON=SSRD;VD=snow surf temperature for roads                           ;VB=p0')
   PHYVAR2D1(sroad_tsen,   'VN=sroad_tsen   ;ON=2MEN; VD=Snow surface temp for roads (E)                          ;VB=e1 ;IN=SSRD;')
   PHYVAR2D1(sroad_wsnow,  'VN=sroad_wsnow  ;ON=SWRD;VD=water and snow content for roads                          ;VB=p0')
   PHYVAR2D1(sroad_wsnowen,'VN=sroad_wsnowen;ON=6ZEN; VD=water/snow content for roads (E)                         ;VB=e1 ;IN=SWRD;')
   PHYVAR2D1(sroof_alb,    'VN=sroof_alb    ;ON=SARF;VD=snow albedo for roofs                                     ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_alben,  'VN=sroof_alben  ;ON=1ZEN; VD=snow albedo for roofs (E)                                ;VB=e1 ;IN=SARF;MIN=0')
   PHYVAR2D1(sroof_emis,   'VN=sroof_emis   ;ON=SERF;VD=snow emissivity for roofs                                 ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_emisen, 'VN=sroof_emisen ;ON=2ZEN; VD=snow emmissivity for roofs (E)                           ;VB=e1 ;IN=SERF;MIN=0')
   PHYVAR2D1(sroof_rho,    'VN=sroof_rho    ;ON=SDRF;VD=snow density for roofs                                    ;VB=p0         ;MIN=0')
   PHYVAR2D1(sroof_rhoen,  'VN=sroof_rhoen  ;ON=9YEN; VD=snow density for roofs (E)                               ;VB=e1 ;IN=SDRF;MIN=0')
   PHYVAR2D1(sroof_scheme, 'VN=sroof_scheme ;ON=SCRF;VD=snow scheme for roofs                                     ;VB=p0')
   PHYVAR2D1(sroof_t,      'VN=sroof_t      ;ON=STRF;VD=snow temperature for roofs                                ;VB=p0')
   PHYVAR2D1(sroof_ten,    'VN=sroof_ten    ;ON=8YEN; VD=snow temp for roofs (E)                                  ;VB=e1 ;IN=STRF;')
   PHYVAR2D1(sroof_ts,     'VN=sroof_ts     ;ON=SSRF;VD=snow surf temperature for roofs                           ;VB=p0')
   PHYVAR2D1(sroof_tsen,   'VN=sroof_tsen   ;ON=3ZEN; VD=Snow surface temp for roofs (E)                          ;VB=e1 ;IN=SSRF;')
   PHYVAR2D1(sroof_wsnow,  'VN=sroof_wsnow  ;ON=SWRF;VD=water and snow content for roofs                          ;VB=p0')
   PHYVAR2D1(sroof_wsnowen,'VN=sroof_wsnowen;ON=7YEN; VD=water/snow content for roofs (E)                         ;VB=e1 ;IN=SWRF;')
   PHYVAR2D1(svf_road,     'VN=svf_road     ;ON=SVRD;VD=road sky-view factor                                      ;VB=p0')
   PHYVAR2D1(svf_wall,     'VN=svf_wall     ;ON=SVWL;VD=wall sky-view factor                                      ;VB=p0')
   PHYVAR2D1(t_canyon,     'VN=t_canyon     ;ON=TCAN;VD=air temperature inside the canyon                         ;VB=p0')
   PHYVAR2D1(t_canyonen,   'VN=t_canyonen   ;ON=1XEN; VD=canyon air temp (E)                                      ;VB=e1 ;IN=TCAN;')
   PHYVAR3D1(t_road,       'VN=t_road       ;ON=TLRD;VD=temperatures of road layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_roaden,     'VN=t_roaden     ;ON=4XEN; VD=road temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLRD;')
   PHYVAR3D1(t_roof,       'VN=t_roof       ;ON=TLRF;VD=temperatures of roof layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_roofen,     'VN=t_roofen     ;ON=3XEN; VD=roof temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLRF;')
   PHYVAR3D1(t_wall,       'VN=t_wall       ;ON=TLWL;VD=temperatures of wall layers                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(t_wallen,     'VN=t_wallen     ;ON=5XEN; VD=wall temperatures (E)                         ;VS=A*3  ;VB=e1 ;IN=TLWL;')
   PHYVAR3D1(tc_road,      'VN=tc_road      ;ON=TCRD;VD=road thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_roaden,    'VN=tc_roaden    ;ON=TB22;VD=road thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCRD;')
   PHYVAR3D1(tc_roof,      'VN=tc_roof      ;ON=TCRF;VD=roof thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_roofen,    'VN=tc_roofen    ;ON=TB23;VD=roof thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCRF;')
   PHYVAR3D1(tc_wall,      'VN=tc_wall      ;ON=TCWL;VD=wall thermal condcutivities                    ;VS=A*3  ;VB=p0')
   PHYVAR3D1(tc_wallen,    'VN=tc_wallen    ;ON=TB24;VD=wall thermal condcutivities (E)                ;VS=A*3  ;VB=e1; IN=TCWL;')
   PHYVAR2D1(ti_bld,       'VN=ti_bld       ;ON=TBLD;VD=internal building temperature                             ;VB=p0')
   PHYVAR2D1(ti_blden,     'VN=ti_blden     ;ON=6QEN; VD=bld internal temp (E)                                    ;VB=e1 ;IN=TBLD;')
   PHYVAR2D1(ti_road,      'VN=ti_road      ;ON=TIRD;VD=internal road temperature                                 ;VB=p0')
   PHYVAR2D1(ti_roaden,    'VN=ti_roaden    ;ON=5QEN; VD=road internal temp (E)                                   ;VB=e1 ;IN=TIRD;')
   PHYVAR2D1(tsradtw,      'VN=tsradtw      ;ON=TSTW;VD=Town radiative surface temperature                        ;VB=p0')
   PHYVAR2D1(tsun,         'VN=tsun         ;ON=TSUN;VD=solar time (s)                                            ;VB=p0')
   PHYVAR2D1(u_canyon,     'VN=u_canyon     ;ON=UCAN;VD=wind in canyon                                            ;VB=p0')
   PHYVAR2D1(wall_o_hor,   'VN=wall_o_hor   ;ON=WHOR;VD=ratio vertical per horizontal surf                        ;VB=p0         ;MIN=0')
   PHYVAR2D1(wall_o_horen, 'VN=wall_o_horen ;ON=TB5 ;VD=ratio vertical per horizontal surf (E)                    ;VB=e1; IN=WHOR;MIN=0')
   PHYVAR2D1(ws_road,      'VN=ws_road      ;ON=WSRD;VD=water content of road reservoir                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(ws_roaden,    'VN=ws_roaden    ;ON=4QEN; VD=road water reservoir (E)                                 ;VB=e1 ;IN=WSRD;MIN=0')
   PHYVAR2D1(ws_roof,      'VN=ws_roof      ;ON=WSRF;VD=water content of roof reservoir                           ;VB=p0         ;MIN=0')
   PHYVAR2D1(ws_roofen,    'VN=ws_roofen    ;ON=3QEN; VD=roof water reservoir (E)                                 ;VB=e1 ;IN=WSRF;MIN=0')

   PHYVAR2DC(yradin,       'VN=yradin       ;ON=RTIN;VD=MRT inside building  (K)                                  ;VB=v0', thermal_stress)
   PHYVAR2DC(yradrfsun,    'VN=yradrfsun    ;ON=RTFS;VD=MRT on the exposed sunny roof (K)                         ;VB=v0', thermal_stress)
   PHYVAR2DC(yradrfshade,  'VN=yradrfshade  ;ON=RTFD;VD=MRT on the shaded roof (K)                                ;VB=v0', thermal_stress)
   PHYVAR2DC(yutciin,      'VN=yutciin      ;ON=DXIN;VD=UTCI inside building  (C)                                 ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcirfsun,   'VN=yutcirfsun   ;ON=DXFS;VD=UTCI on the exposed sunny roof (C)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcirfshade, 'VN=yutcirfshade ;ON=DXFD;VD=UTCI on the shaded roof (C)                               ;VB=v0', thermal_stress)
   PHYVAR2DC(ytrfzt,       'VN=ytrfzt       ;ON=T2RF;VD=Temperature at zt above the roof                          ;VB=v0', thermal_stress)
   PHYVAR2DC(ytrdzt,       'VN=ytrdzt       ;ON=T2RD;VD=Temperature at zt above the ground                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yurdzu,       'VN=yurdzu       ;ON=UVRD;VD=wind speed at zu above the ground (m/s)                   ;VB=v0', thermal_stress)
   PHYVAR2DC(ywbgtrfsun,   'VN=ywbgtrfsun   ;ON=GXFS;VD=WBGT on the exposed sunny roof (C)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(ywbgtrfshade, 'VN=ywbgtrfshade ;ON=GXFD;VD=WBGT on the shaded roof (C)                               ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicin,     'VN=yutcicin     ;ON=DCIN;VD=cumulative UTCI inside building  (C)                      ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicsun,    'VN=yutcicsun    ;ON=DCSU;VD=cumulative UTCI in the exposed sunny street               ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicshade,  'VN=yutcicshade  ;ON=DCHD;VD=cumulative UTCI in the shaded street (C)                  ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicrfsun,  'VN=yutcicrfsun  ;ON=DCFS;VD=cumulative UTCI on the exposed sunny roof (C)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yutcicrfshade,'VN=yutcicrfshade;ON=DCFD;VD=cumulative UTCI on the shaded roof (C)                    ;VB=v0', thermal_stress)
   PHYVAR2DC(ytglbrfsun,   'VN=ytglbrfsun   ;ON=GTFS;VD=TGlobe on the exposed sunny roof (K)                      ;VB=v0', thermal_stress)
   PHYVAR2DC(ytglbrfshade, 'VN=ytglbrfshade ;ON=GTFD;VD=TGlobe on the shaded roof (K)                             ;VB=v0', thermal_stress)
   PHYVAR2DC(ytwetbrf,     'VN=ytwetbrf     ;ON=WBRF;VD=Wet-Bulb Temperature at zt above the roof (K)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ8,          'VN=yQ8          ;ON=QSRF;VD=Contribution of roof SW rad (W/m2)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ9,          'VN=yQ9          ;ON=QLRF;VD=Contribution of roof LW rad (W/m2)                        ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ10,         'VN=yQ10         ;ON=QSFK;VD=Contribution of sky on the roof SW rad (W/m2)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ11,         'VN=yQ11         ;ON=QLFK;VD=Contribution of sky on the roof LW rad (W/m2)             ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ12,         'VN=yQ12         ;ON=QSFW;VD=Contribution of wall on the roof SW rad (W/m2)            ;VB=v0', thermal_stress)
   PHYVAR2DC(yQ13,         'VN=yQ13         ;ON=QLFW;VD=Contribution of wall on the roof LW rad (W/m2)            ;VB=v0', thermal_stress)

   PHYVAR2D1(z0_road,      'VN=z0_road      ;ON=Z0RD;VD=aerodyn roughness length for road                         ;VB=p0')
   PHYVAR2D1(z0_roaden,    'VN=z0_roaden    ;ON=TB8 ;VD=aerodyn roughness length for road (E)                     ;VB=e1; IN=Z0RD;')
   PHYVAR2D1(z0_roof,      'VN=z0_roof      ;ON=Z0RF;VD=aerodyn roughness length for roof                         ;VB=p0')
   PHYVAR2D1(z0_roofen,    'VN=z0_roofen    ;ON=TB7 ;VD=aerodyn roughness length for roof (E)                     ;VB=e1; IN=Z0RF;')
   PHYVAR2D1(z0_town,      'VN=z0_town      ;ON=Z0TW;VD=aerodyn roughness length for town                         ;VB=p0')
   PHYVAR2D1(z0_townen,    'VN=z0_townen    ;ON=TB6 ;VD=aerodyn roughness length for town (E)                     ;VB=e1; IN=Z0TW;')
   PHYVAR2D1(zenith,       'VN=zenith       ;ON=ZENI;VD=solar zenith angle                                        ;VB=p0')

   !---------------------------------------------------------------------
   return
end subroutine sfc_businit
