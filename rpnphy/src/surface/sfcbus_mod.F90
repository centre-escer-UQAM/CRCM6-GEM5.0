!-------------------------------------- LICENCE BEGIN ------------------------
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
!-------------------------------------- LICENCE END --------------------------

module sfcbus_mod
   implicit none
   public
   save

#define SFCVAR(MYVAR,MYNAME) type(SFCVAR_T) :: MYVAR = SFCVAR_T(-1,0,0,0,0,.false.,MYNAME)

   type :: sfcptr
      sequence
      real, pointer :: ptr(:,:)
   end type sfcptr

   type :: SFCVAR_T
      sequence
      integer :: i, agg, mul, niveaux, mosaik
      logical :: doagg_L
      character(len=32) :: n
   end type SFCVAR_T

   type :: SFCVARLIST_T
      sequence

      SFCVAR(umoins , 'PW_UU:M')
      SFCVAR(uplus  , 'PW_UU:P')
      SFCVAR(vmoins , 'PW_VV:M')
      SFCVAR(vplus  , 'PW_VV:P')
      SFCVAR(tmoins , 'PW_TT:M')
      SFCVAR(tplus  , 'PW_TT:P')
      SFCVAR(sigm   , 'PW_PM:P')
      SFCVAR(humoins, 'TR/HU:M')
      SFCVAR(huplus , 'TR/HU:P')

      SFCVAR(acoef, 'acoef')
      SFCVAR(acroot, 'acroot')
      SFCVAR(accevap, 'accevap')
      SFCVAR(alb_road, 'alb_road')
      SFCVAR(alb_roaden, 'alb_roaden')
      SFCVAR(alb_roof, 'alb_roof')
      SFCVAR(alb_roofen, 'alb_roofen')
      SFCVAR(alb_wall, 'alb_wall')
      SFCVAR(alb_wallen, 'alb_wallen')
      SFCVAR(alfaq, 'alfaq')
      SFCVAR(alfat, 'alfat')
      SFCVAR(algr, 'algr')
      SFCVAR(alir, 'alir')
      SFCVAR(alscatw, 'alscatw')
      SFCVAR(alveg, 'alveg')
      SFCVAR(alvh, 'alvh')
      SFCVAR(alvis, 'alvis')
      SFCVAR(alvl, 'alvl')
      SFCVAR(alvs, 'alvs')
      SFCVAR(algdry, 'algdry')
      SFCVAR(algwet, 'algwet')
      SFCVAR(alwater, 'alwater')
      SFCVAR(avg_gwsol, 'avg_gwsol')
      SFCVAR(azim, 'azim')
      SFCVAR(bcoef, 'bcoef')
      SFCVAR(bld, 'bld')
      SFCVAR(blden, 'blden')
      SFCVAR(bld_height, 'bld_height')
      SFCVAR(bld_heighten, 'bld_heighten')
      SFCVAR(bm, 'bm')
      SFCVAR(bt, 'bt')
      SFCVAR(c1sat, 'c1sat')
      SFCVAR(c2ref, 'c2ref')
      SFCVAR(c3ref, 'c3ref')
      SFCVAR(cang, 'cang')
      SFCVAR(can_hw_ratio, 'can_hw_ratio')
      SFCVAR(cgsat, 'cgsat')
      SFCVAR(clay, 'clay')
      SFCVAR(clayen, 'clayen')
      SFCVAR(co2i1, 'co2i1')
      SFCVAR(cosz, 'cosz')
      SFCVAR(cveg, 'cveg')
      SFCVAR(cvh, 'cvh')
      SFCVAR(cvl, 'cvl')
      SFCVAR(d50, 'd50')
      SFCVAR(d95, 'd95')
      SFCVAR(d_road, 'd_road')
      SFCVAR(d_roaden, 'd_roaden')
      SFCVAR(d_roof, 'd_roof')
      SFCVAR(d_roofen, 'd_roofen')
      SFCVAR(d_wall, 'd_wall')
      SFCVAR(d_wallen, 'd_wallen')
      SFCVAR(deciduous, 'deciduous')
      SFCVAR(dhdx, 'dhdx')
      SFCVAR(dhdxdy, 'dhdxdy')
      SFCVAR(dhdy, 'dhdy')
      SFCVAR(dlat, 'dlat')
      SFCVAR(dlon, 'dlon')
      SFCVAR(drain, 'drain')
      SFCVAR(draindens, 'draindens')
      SFCVAR(drainaf, 'drainaf')
      SFCVAR(draintot, 'draintot')
      SFCVAR(draintotaf, 'draintotaf')
      SFCVAR(dsst, 'dsst')
      SFCVAR(dtdiag, 'dtdiag')
      SFCVAR(eflux, 'eflux')
      SFCVAR(eg, 'eg')
      SFCVAR(emis, 'emis')
      SFCVAR(emisr, 'emisr')
      SFCVAR(emis_road, 'emis_road')
      SFCVAR(emis_roaden, 'emis_roaden')
      SFCVAR(emis_roof, 'emis_roof')
      SFCVAR(emis_roofen, 'emis_roofen')
      SFCVAR(emis_wall, 'emis_wall')
      SFCVAR(emis_wallen, 'emis_wallen')
      SFCVAR(emisgr, 'emisgr')
      SFCVAR(emistg, 'emistg')
      SFCVAR(emistgen, 'emistgen')
      SFCVAR(emsvc, 'emsvc')
      SFCVAR(emisvh, 'emisvh')
      SFCVAR(emisvl, 'emisvl')
      SFCVAR(emtw, 'emtw')
      SFCVAR(en, 'en')
      SFCVAR(er, 'er')
      SFCVAR(etr, 'etr')
      SFCVAR(evergreen, 'evergreen')
      SFCVAR(fbcof, 'fbcof')
      SFCVAR(fc, 'fc')
      SFCVAR(fcor, 'fcor')
      SFCVAR(fdsi, 'fdsi')
      SFCVAR(fdss, 'fdss')
      SFCVAR(fl, 'fl')
      SFCVAR(flusolis, 'flusolis')
      SFCVAR(fluslop, 'fluslop')
      SFCVAR(fq, 'fq')
      SFCVAR(frootd, 'frootd')
      SFCVAR(frv, 'frv')
      SFCVAR(frv_li, 'frv_li')
      SFCVAR(frv_lw, 'frv_lw')
      SFCVAR(fsd, 'fsd')
      SFCVAR(fsf, 'fsf')
      SFCVAR(ftemp, 'ftemp')
      SFCVAR(fv, 'fv')
      SFCVAR(fvap, 'fvap')
      SFCVAR(fvapliq, 'fvapliq')
      SFCVAR(fvapliqaf, 'fvapliqaf')
      SFCVAR(g_road, 'g_road')
      SFCVAR(g_roof, 'g_roof')
      SFCVAR(g_town, 'g_town')
      SFCVAR(g_wall, 'g_wall')
      SFCVAR(gamveg, 'gamveg')
      SFCVAR(gamvh, 'gamvh')
      SFCVAR(gamvl, 'gamvl')      
      SFCVAR(gc, 'gc')
      SFCVAR(glacier, 'glacier')
      SFCVAR(glsea, 'glsea')
      SFCVAR(glsea0, 'glsea0')
      SFCVAR(grkef, 'grkef')
      SFCVAR(grksat, 'grksat')
      SFCVAR(gztherm, 'gztherm')
      SFCVAR(h, 'h')
      SFCVAR(h_industry, 'h_industry')
      SFCVAR(h_industryen, 'h_industryen')
      SFCVAR(h_road, 'h_road')
      SFCVAR(h_roof, 'h_roof')
      SFCVAR(h_town, 'h_town')
      SFCVAR(h_traffic, 'h_traffic')
      SFCVAR(h_trafficen, 'h_trafficen')
      SFCVAR(h_wall, 'h_wall')
      SFCVAR(hc_road, 'hc_road')
      SFCVAR(hc_roaden, 'hc_roaden')
      SFCVAR(hc_roof, 'hc_roof')
      SFCVAR(hc_roofen, 'hc_roofen')
      SFCVAR(hc_wall, 'hc_wall')
      SFCVAR(hc_wallen, 'hc_wallen')
      SFCVAR(hfluxsa, 'hfluxsa')
      SFCVAR(hfluxsv, 'hfluxsv')
      SFCVAR(hst, 'hst')
      SFCVAR(husurf, 'husurf')
      SFCVAR(hv, 'hv')
      SFCVAR(icedp, 'icedp')
      SFCVAR(iceline, 'iceline')
      SFCVAR(ilmo, 'ilmo')
      SFCVAR(impervu, 'impervu')
      SFCVAR(isoil, 'isoil')
      SFCVAR(kcl, 'kcl')
      SFCVAR(khc, 'khc')
      SFCVAR(km, 'km')
      SFCVAR(ksat, 'ksat')
      SFCVAR(ksatc, 'ksatc')
      SFCVAR(kt, 'kt')
      SFCVAR(lai, 'lai')
      SFCVAR(laictem, 'laictem')
      SFCVAR(laideci, 'laideci')
      SFCVAR(laiva, 'laiva')
      SFCVAR(laivf26, 'laivf26')
      SFCVAR(laivh, 'laivh')
      SFCVAR(laivl, 'laivl')
      SFCVAR(lakect, 'lakect')
      SFCVAR(lakedepth, 'lakedepth')
      SFCVAR(lakefice, 'lakefice')
      SFCVAR(lakefr, 'lakefr')
      SFCVAR(lakehice, 'lakehice')
      SFCVAR(lakehml, 'lakehml')
      SFCVAR(laketbot, 'laketbot')
      SFCVAR(laketice, 'laketice')
      SFCVAR(laketmnw, 'laketmnw')
      SFCVAR(laketp, 'laketp')
      SFCVAR(laketwml, 'laketwml')
      SFCVAR(laketransp, 'laketransp')
      SFCVAR(latflaf, 'latflaf')
      SFCVAR(latflw, 'latflw')
      SFCVAR(le_industry, 'le_industry')
      SFCVAR(le_industryen, 'le_industryen')
      SFCVAR(le_road, 'le_road')
      SFCVAR(le_roof, 'le_roof')
      SFCVAR(le_town, 'le_town')
      SFCVAR(le_traffic, 'le_traffic')
      SFCVAR(le_trafficen, 'le_trafficen')
      SFCVAR(le_wall, 'le_wall')
      SFCVAR(leg, 'leg')
      SFCVAR(ler, 'ler')
      SFCVAR(les, 'les')
      SFCVAR(lesv, 'lesv')
      SFCVAR(letr, 'letr')
      SFCVAR(lev, 'lev')
      SFCVAR(lhtg, 'lhtg')
      SFCVAR(melts, 'melts')
      SFCVAR(meltsr, 'meltsr')
      SFCVAR(mf, 'mf')
      SFCVAR(mg, 'mg')
      SFCVAR(ml, 'ml')
      SFCVAR(mt, 'mt')
      SFCVAR(mtdir, 'mtdir')
      SFCVAR(nat, 'nat')
      SFCVAR(overfl, 'overfl')
      SFCVAR(pav, 'pav')
      SFCVAR(paven, 'paven')
      SFCVAR(pcoef, 'pcoef')
      SFCVAR(pmoins, 'pmoins')
      SFCVAR(pplus, 'pplus')
      SFCVAR(psi, 'psi')
      SFCVAR(psisat, 'psisat')
      SFCVAR(psn, 'psn')
      SFCVAR(psng, 'psng')
      SFCVAR(psngrvl, 'psngrvl')
      SFCVAR(psnv, 'psnv')
      SFCVAR(psnvh, 'psnvh')
      SFCVAR(psnvha, 'psnvha')
      SFCVAR(q_canyon, 'q_canyon')
      SFCVAR(q_canyonen, 'q_canyonen')
      SFCVAR(qdiag, 'qdiag')
      SFCVAR(qdiagstn, 'qdiagstn')
      SFCVAR(qdiagstnv, 'qdiagstnv')
      SFCVAR(qdiagtyp, 'qdiagtyp')
      SFCVAR(qdiagtypv, 'qdiagtypv')
      SFCVAR(qsurf, 'qsurf')
      SFCVAR(rainrate, 'rainrate')
      SFCVAR(rcctem, 'rcctem')
      SFCVAR(resa, 'resa')
      SFCVAR(resagr, 'resagr')
      SFCVAR(resavg, 'resavg')
      SFCVAR(resasa, 'resasa')
      SFCVAR(resasv, 'resasv')
      SFCVAR(resaef, 'resaef')
      SFCVAR(rgl, 'rgl')
      SFCVAR(rglvh, 'rglvh')
      SFCVAR(rglvl, 'rglvl')
      SFCVAR(rib, 'rib')
      SFCVAR(rn_road, 'rn_road')
      SFCVAR(rn_roof, 'rn_roof')
      SFCVAR(rn_town, 'rn_town')
      SFCVAR(rn_wall, 'rn_wall')
      SFCVAR(rnet_s, 'rnet_s')
      SFCVAR(rnetsa, 'rnetsa')
      SFCVAR(rnetsv, 'rnetsv')
      SFCVAR(rootdp, 'rootdp')
      SFCVAR(rsnowsa, 'rsnowsa')
      SFCVAR(rsnowsv, 'rsnowsv')
      SFCVAR(rst, 'rst')
      SFCVAR(rt, 'rt')
      SFCVAR(runofftot, 'runofftot')
      SFCVAR(runofftotaf, 'runofftotaf')
      SFCVAR(rveg, 'rveg')
      SFCVAR(sand, 'sand')
      SFCVAR(sanden, 'sanden')
      SFCVAR(sfcwgt, 'sfcwgt')
      SFCVAR(skin_depth, 'skin_depth')
      SFCVAR(skin_inc, 'skin_inc')
      SFCVAR(skyview, 'skyview')
      SFCVAR(slop, 'slop')
      SFCVAR(slope, 'slope')
      SFCVAR(snden, 'snden')
      SFCVAR(snoagen, 'snoagen')
      SFCVAR(snoal, 'snoal')
      SFCVAR(snoalen, 'snoalen')
      SFCVAR(snoden, 'snoden')
      SFCVAR(snodp, 'snodp')
      SFCVAR(snodpl, 'snodpl')
      SFCVAR(snoma, 'snoma')
      SFCVAR(snoro, 'snoro')
      SFCVAR(snowrate, 'snowrate')
      SFCVAR(snowsize, 'snowsize')
      SFCVAR(snval, 'snval')
      SFCVAR(snvden, 'snvden')
      SFCVAR(snvdp, 'snvdp')
      SFCVAR(snvma, 'snvma')
      SFCVAR(snvro, 'snvro')  
      SFCVAR(sroad_alb, 'sroad_alb')
      SFCVAR(sroad_alben, 'sroad_alben')
      SFCVAR(sroad_emis, 'sroad_emis')
      SFCVAR(sroad_emisen, 'sroad_emisen')
      SFCVAR(sroad_rho, 'sroad_rho')
      SFCVAR(sroad_rhoen, 'sroad_rhoen')
      SFCVAR(sroad_t, 'sroad_t')
      SFCVAR(sroad_ten, 'sroad_ten')
      SFCVAR(sroad_ts, 'sroad_ts')
      SFCVAR(sroad_tsen, 'sroad_tsen')
      SFCVAR(sroad_wsnow, 'sroad_wsnow')
      SFCVAR(sroad_wsnowen, 'sroad_wsnowen')
      SFCVAR(sroof_alb, 'sroof_alb')
      SFCVAR(sroof_alben, 'sroof_alben')
      SFCVAR(sroof_emis, 'sroof_emis')
      SFCVAR(sroof_emisen, 'sroof_emisen')
      SFCVAR(sroof_rho, 'sroof_rho')
      SFCVAR(sroof_rhoen, 'sroof_rhoen')
      SFCVAR(sroof_t, 'sroof_t')
      SFCVAR(sroof_ten, 'sroof_ten')
      SFCVAR(sroof_ts, 'sroof_ts')
      SFCVAR(sroof_tsen, 'sroof_tsen')
      SFCVAR(sroof_wsnow, 'sroof_wsnow')
      SFCVAR(sroof_wsnowen, 'sroof_wsnowen')
      SFCVAR(stomr, 'stomr')
      SFCVAR(stomrvh, 'stomrvh')
      SFCVAR(stomrvl, 'stomrvl')
      SFCVAR(svf_road, 'svf_road')
      SFCVAR(svf_wall, 'svf_wall')
      SFCVAR(svs_wta, 'svs_wta')    
      SFCVAR(sw4diff, 'sw4diff')    
      SFCVAR(sw4drct, 'sw4drct')    
      SFCVAR(sw4totl, 'sw4totl')    
      SFCVAR(t_canyon, 't_canyon')
      SFCVAR(t_canyonen, 't_canyonen')
      SFCVAR(t_road, 't_road')
      SFCVAR(t_roaden, 't_roaden')
      SFCVAR(t_roof, 't_roof')
      SFCVAR(t_roofen, 't_roofen')
      SFCVAR(t_wall, 't_wall')
      SFCVAR(t_wallen, 't_wallen')
      SFCVAR(tc_road, 'tc_road')
      SFCVAR(tc_roaden, 'tc_roaden')
      SFCVAR(tc_roof, 'tc_roof')
      SFCVAR(tc_roofen, 'tc_roofen')
      SFCVAR(tc_wall, 'tc_wall')
      SFCVAR(tc_wallen, 'tc_wallen')
      SFCVAR(tddiagtyp, 'tddiagtyp')
      SFCVAR(tddiagtypv, 'tddiagtypv')
      SFCVAR(tdiag, 'tdiag')
      SFCVAR(tdiagstn, 'tdiagstn')
      SFCVAR(tdiagstnv, 'tdiagstnv')
      SFCVAR(tdiagtyp, 'tdiagtyp')
      SFCVAR(tdiagtypv, 'tdiagtypv')
      SFCVAR(tglacier, 'tglacier')
      SFCVAR(tground, 'tground')
      SFCVAR(thetaa, 'thetaa')
      SFCVAR(thetaap, 'thetaap')
      SFCVAR(ti_bld, 'ti_bld')
      SFCVAR(ti_blden, 'ti_blden')
      SFCVAR(ti_road, 'ti_road')
      SFCVAR(ti_roaden, 'ti_roaden')
      SFCVAR(tmice, 'tmice')
      SFCVAR(tnolim, 'tnolim')
      SFCVAR(trunofftot, 'trunofftot')
      SFCVAR(trunofftotaf, 'trunofftotaf')
      SFCVAR(tsa, 'tsa')
      SFCVAR(tsnavg, 'tsnavg')
      SFCVAR(tsnow, 'tsnow')
      SFCVAR(tsnowveg, 'tsnowveg')
      SFCVAR(tsoil, 'tsoil')
      SFCVAR(tsrad, 'tsrad')
      SFCVAR(tsradtw, 'tsradtw')
      SFCVAR(tss, 'tss')
      SFCVAR(tsun, 'tsun')
      SFCVAR(tsurf, 'tsurf')
      SFCVAR(tsvavg, 'tsvavg')
      SFCVAR(tve, 'tve')
      SFCVAR(tvege, 'tvege')
      SFCVAR(twater, 'twater')
      SFCVAR(u_canyon, 'u_canyon')
      SFCVAR(udiag, 'udiag')
      SFCVAR(udiagstn, 'udiagstn')
      SFCVAR(udiagstnv, 'udiagstnv')
      SFCVAR(udiagtyp, 'udiagtyp')
      SFCVAR(udiagtypv, 'udiagtypv')
      SFCVAR(urban, 'urban')
      SFCVAR(vdiag, 'vdiag')
      SFCVAR(vdiagstn, 'vdiagstn')
      SFCVAR(vdiagstnv, 'vdiagstnv')
      SFCVAR(vdiagtyp, 'vdiagtyp')
      SFCVAR(vdiagtypv, 'vdiagtypv')
      SFCVAR(vegf, 'vegf')
      SFCVAR(vegfrac, 'vegfrac')
      SFCVAR(vegh, 'vegh')
      SFCVAR(vegl, 'vegl')
      SFCVAR(vegtrans, 'vegtrans')
      SFCVAR(vgctem, 'vgctem')
      SFCVAR(wall_o_hor, 'wall_o_hor')
      SFCVAR(wall_o_horen, 'wall_o_horen')
      SFCVAR(watflow, 'watflow')
      SFCVAR(wfc, 'wfc')
      SFCVAR(wfcdp, 'wfcdp')
      SFCVAR(wfcint, 'wfcint')
      SFCVAR(wflux, 'wflux')
      SFCVAR(ws_road, 'ws_road')
      SFCVAR(ws_roaden, 'ws_roaden')
      SFCVAR(ws_roof, 'ws_roof')
      SFCVAR(ws_roofen, 'ws_roofen')
      SFCVAR(wsat, 'wsat')
      SFCVAR(wsnow, 'wsnow')
      SFCVAR(wsnv, 'wsnv')
      SFCVAR(wsoil, 'wsoil')
      SFCVAR(wsoilm, 'wsoilm')
      SFCVAR(wveg, 'wveg')
      SFCVAR(wwilt, 'wwilt')
      SFCVAR(xcent, 'xcent')
      SFCVAR(yradin, 'yradin')
      SFCVAR(yradsun, 'yradsun')
      SFCVAR(yradshade, 'yradshade')
      SFCVAR(yradrfsun, 'yradrfsun')
      SFCVAR(yradrfshade, 'yradrfshade')
      SFCVAR(yutciin, 'yutciin')
      SFCVAR(yutcisun, 'yutcisun')
      SFCVAR(yutcishade, 'yutcishade')
      SFCVAR(yutcirfsun, 'yutcirfsun')
      SFCVAR(yutcirfshade, 'yutcirfshade')
      SFCVAR(ywbgtsun, 'ywbgtsun')
      SFCVAR(ywbgtshade, 'ywbgtshade')
      SFCVAR(ywbgtrfsun, 'ywbgtrfsun')
      SFCVAR(ywbgtrfshade, 'ywbgtrfshade')
      SFCVAR(yutcicin, 'yutcicin')
      SFCVAR(yutcicsun, 'yutcicsun')
      SFCVAR(yutcicshade, 'yutcicshade')
      SFCVAR(yutcicrfsun, 'yutcicrfsun')
      SFCVAR(yutcicrfshade, 'yutcicrfshade')
      SFCVAR(ytglbsun, 'ytglbsun')
      SFCVAR(ytglbshade, 'ytglbshade')
      SFCVAR(ytglbrfsun, 'ytglbrfsun')
      SFCVAR(ytglbrfshade, 'ytglbrfshade')
      SFCVAR(ytwetb, 'ytwetb')
      SFCVAR(ytwetbrf, 'ytwetbrf')
      SFCVAR(ytrfzt, 'ytrfzt')
      SFCVAR(ytrdzt, 'ytrdzt')
      SFCVAR(yurdzu, 'yurdzu')
      SFCVAR(yQ1, 'yQ1')
      SFCVAR(yQ2, 'yQ2')
      SFCVAR(yQ3, 'yQ3')
      SFCVAR(yQ4, 'yQ4')
      SFCVAR(yQ5, 'yQ5')
      SFCVAR(yQ6, 'yQ6')
      SFCVAR(yQ7, 'yQ7')
      SFCVAR(yQ8, 'yQ8')
      SFCVAR(yQ9, 'yQ9')
      SFCVAR(yQ10, 'yQ10')
      SFCVAR(yQ11, 'yQ11')
      SFCVAR(yQ12, 'yQ12')
      SFCVAR(yQ13, 'yQ13')
      SFCVAR(z0, 'z0')
      SFCVAR(z0_road, 'z0_road')
      SFCVAR(z0_roaden, 'z0_roaden')
      SFCVAR(z0_roof, 'z0_roof')
      SFCVAR(z0_roofen, 'z0_roofen')
      SFCVAR(z0_town, 'z0_town')
      SFCVAR(z0_townen, 'z0_townen')
      SFCVAR(z0en, 'z0en')
      SFCVAR(z0ha, 'z0ha')
      SFCVAR(z0mvh, 'z0mvh')
      SFCVAR(z0mvl, 'z0mvl')
      SFCVAR(z0veg, 'z0veg')
      SFCVAR(z0t, 'z0t')
      SFCVAR(z0tveg, 'z0tveg')
      SFCVAR(za, 'za')
      SFCVAR(ze, 'ze')
      SFCVAR(zenith, 'zenith')
      SFCVAR(ztsl, 'ztsl')
      SFCVAR(zusl, 'zusl')

      ! CLASS fields
      SFCVAR(algdn, 'algdn')
      SFCVAR(algdv, 'algdv')
      SFCVAR(algwn, 'algwn')
      SFCVAR(algwv, 'algwv')
      SFCVAR(alirc, 'alirc')
      SFCVAR(alvsc, 'alvsc')
      SFCVAR(bbi, 'bbi')
      SFCVAR(cdh, 'cdh')
      SFCVAR(cdm, 'cdm')
      SFCVAR(cmai, 'cmai')
      SFCVAR(delzw, 'delzw')
      SFCVAR(evapo, 'evapo')
      SFCVAR(fcanmx, 'fcanmx')
      SFCVAR(fcovc, 'fcovc')
      SFCVAR(fcovcs, 'fcovcs')
      SFCVAR(fcovg, 'fcovg')
      SFCVAR(fcovgs, 'fcovgs')
      SFCVAR(firupaf, 'firupaf')
      SFCVAR(flgg, 'flgg')
      SFCVAR(flgs, 'flgs')
      SFCVAR(flgv, 'flgv')
      SFCVAR(fsgg, 'fsgg')
      SFCVAR(fsgs, 'fsgs')
      SFCVAR(fsgv, 'fsgv')
      SFCVAR(fsnow, 'fsnow')
      SFCVAR(fsolupaf, 'fsolupaf')
      SFCVAR(grdhflx, 'grdhflx')
      SFCVAR(grkfac, 'grkfac')
      SFCVAR(hcps, 'hcps')
      SFCVAR(hevc, 'hevc')
      SFCVAR(hevg, 'hevg')
      SFCVAR(hevs, 'hevs')
      SFCVAR(hfsc, 'hfsc')
      SFCVAR(hfsg, 'hfsg')
      SFCVAR(hfss, 'hfss')
      SFCVAR(hmfc, 'hmfc')
      SFCVAR(hmfg, 'hmfg')
      SFCVAR(hmfn, 'hmfn')
      SFCVAR(htc, 'htc')
      SFCVAR(htcc, 'htcc')
      SFCVAR(htcs, 'htcs')
      SFCVAR(huaircan, 'huaircan')
      SFCVAR(iveg, 'iveg')
      SFCVAR(laimax, 'laimax')
      SFCVAR(laimin, 'laimin')
      SFCVAR(mosfract, 'mosfract')
      SFCVAR(orgm, 'orgm')
      SFCVAR(pcfc, 'pcfc')
      SFCVAR(pcpn, 'pcpn')
      SFCVAR(pclc, 'pclc')
      SFCVAR(pcpg, 'pcpg')
      SFCVAR(psiga, 'psiga')
      SFCVAR(psigb, 'psigb')
      SFCVAR(psiwlt, 'psiwlt')
      SFCVAR(qa50, 'qa50')
      SFCVAR(qfc, 'qfc')
      SFCVAR(qfcf, 'qfcf')
      SFCVAR(qfcl, 'qfcl')
      SFCVAR(qfg, 'qfg')
      SFCVAR(qfn, 'qfn')
      SFCVAR(rofc, 'rofc')
      SFCVAR(rofn, 'rofn')
      SFCVAR(rovg, 'rovg')
      SFCVAR(sdepth, 'sdepth')
      SFCVAR(subflw, 'subflw')
      SFCVAR(taircan, 'taircan')
      SFCVAR(tbase, 'tbase')
      SFCVAR(tbasfl, 'tbasfl')
      SFCVAR(tcs, 'tcs')
      SFCVAR(thfc, 'thfc')
      SFCVAR(thlmin, 'thlmin')
      SFCVAR(thlrat, 'thlrat')
      SFCVAR(thlret, 'thlret')
      SFCVAR(thpor, 'thpor')
      SFCVAR(tovrfl, 'tovrfl')
      SFCVAR(tpond, 'tpond')
      SFCVAR(trunoff, 'trunoff')
      SFCVAR(tsno, 'tsno')
      SFCVAR(tsubfl, 'tsubfl')
      SFCVAR(tsurfsa, 'tsurfsa')
      SFCVAR(tveg, 'tveg')
      SFCVAR(veggro, 'veggro')
      SFCVAR(vegma, 'vegma')
      SFCVAR(vpda, 'vpda')
      SFCVAR(vpdb, 'vpdb')
      SFCVAR(wfsurf, 'wfsurf')
      SFCVAR(wtrc, 'wtrc')
      SFCVAR(wtrg, 'wtrg')
      SFCVAR(wtrs, 'wtrs')
      SFCVAR(xdrain, 'xdrain')
      SFCVAR(xslope, 'xslope')
      SFCVAR(z0oro, 'z0oro')
      SFCVAR(z0vegc, 'z0vegc')
      SFCVAR(zbotw, 'zbotw')
      SFCVAR(zoln, 'zoln')
      SFCVAR(zpond, 'zpond')
      SFCVAR(anis, 'anis')
      SFCVAR(are, 'are')
      SFCVAR(excw, 'excw')
      SFCVAR(lbedr, 'lbedr')
      SFCVAR(leggw, 'leggw')
      SFCVAR(slpgw, 'slpgw')
      SFCVAR(totw, 'totw')
      SFCVAR(wtnew, 'wtnew')
      ! End of CLASS fields

      ! CTEM fields
      SFCVAR(ailc, 'ailc')
      SFCVAR(ailcb, 'ailcb')
      SFCVAR(ailcg, 'ailcg')
      SFCVAR(alirctm, 'alirctm')
      SFCVAR(allwacc, 'allwacc')
      SFCVAR(alswacc, 'alswacc')
      SFCVAR(alvsctm, 'alvsctm')
      SFCVAR(ancgvgac, 'ancgvgac')
      SFCVAR(ancsvgac, 'ancsvgac')
      SFCVAR(anndefct, 'anndefct')
      SFCVAR(annpcp, 'annpcp')
      SFCVAR(annsrpls, 'annsrpls')
      SFCVAR(anpcpcur, 'anpcpcur')
      SFCVAR(anpecur, 'anpecur')
      SFCVAR(anpotevp, 'anpotevp')
      SFCVAR(aridity, 'aridity')
      SFCVAR(bleafmas, 'bleafmas')
      SFCVAR(bmasveg, 'bmasveg')
      SFCVAR(burnvegf, 'burnvegf')
      SFCVAR(cfluxcg, 'cfluxcg')
      SFCVAR(cfluxcs, 'cfluxcs')
      SFCVAR(cmasvegc, 'cmasvegc')
      SFCVAR(co2conc, 'co2conc')
      SFCVAR(co2i1cg, 'co2i1cg')
      SFCVAR(co2i1cs, 'co2i1cs')
      SFCVAR(co2i2cg, 'co2i2cg')
      SFCVAR(co2i2cs, 'co2i2cs')
      SFCVAR(colddayr, 'colddayr')
      SFCVAR(defctcur, 'defctcur')
      SFCVAR(defctmon, 'defctmon')
      SFCVAR(defmnr, 'defmnr')
      SFCVAR(dftcuryr, 'dftcuryr')
      SFCVAR(dryslen, 'dryslen')
      SFCVAR(dvdfcan, 'dvdfcan')
      SFCVAR(extnprob, 'extnprob')
      SFCVAR(fcancmx, 'fcancmx')
      SFCVAR(flhrloss, 'flhrloss')
      SFCVAR(flinacc, 'flinacc')
      SFCVAR(flutacc, 'flutacc')
      SFCVAR(fsinacc, 'fsinacc')
      SFCVAR(fsnowacc, 'fsnowacc')
      SFCVAR(gavglai, 'gavglai')
      SFCVAR(gavgltms, 'gavgltms')
      SFCVAR(gavgscms, 'gavgscms')
      SFCVAR(gdd5, 'gdd5')
      SFCVAR(gdd5cur, 'gdd5cur')
      SFCVAR(geremort, 'geremort')
      SFCVAR(gleafmas, 'gleafmas')
      SFCVAR(grwtheff, 'grwtheff')
      SFCVAR(intrmort, 'intrmort')
      SFCVAR(lambda, 'lambda')
      SFCVAR(lfstatur, 'lfstatur')
      SFCVAR(litrmass, 'litrmass')
      SFCVAR(lyrotmas, 'lyrotmas')
      SFCVAR(lystmmas, 'lystmmas')
      SFCVAR(mlightng, 'mlightng')
      SFCVAR(nfcancmx, 'nfcancmx')
      SFCVAR(nppveg, 'nppveg')
      SFCVAR(paic, 'paic')
      SFCVAR(pandayr, 'pandayr')
      SFCVAR(pfcancmx, 'pfcancmx')
      SFCVAR(pftexistr, 'pftexistr')
      SFCVAR(pgleafmass, 'pgleafmass')
      SFCVAR(prbfrhuc, 'prbfrhuc')
      SFCVAR(preacc, 'preacc')
      SFCVAR(pstemmass, 'pstemmass')
      SFCVAR(rmatc, 'rmatc')
      SFCVAR(rmatctem, 'rmatctem')
      SFCVAR(rmlcgvga, 'rmlcgvga')
      SFCVAR(rmlcsvga, 'rmlcsvga')
      SFCVAR(rootdpth, 'rootdpth')
      SFCVAR(rootmass, 'rootmass')
      SFCVAR(rothrlos, 'rothrlos')
      SFCVAR(slai, 'slai')
      SFCVAR(slaic, 'slaic')
      SFCVAR(soilcmas, 'soilcmas')
      SFCVAR(srpcuryr, 'srpcuryr')
      SFCVAR(srplscur, 'srplscur')
      SFCVAR(srplsmon, 'srplsmon')
      SFCVAR(stemmass, 'stemmass')
      SFCVAR(stmhrlos, 'stmhrlos')
      SFCVAR(surmnr, 'surmnr')
      SFCVAR(taaccgat, 'taaccgat')
      SFCVAR(tbaraccgat, 'tbaraccgat')
      SFCVAR(tbarcacc, 'tbarcacc')
      SFCVAR(tbarcsacc, 'tbarcsacc')
      SFCVAR(tbargacc, 'tbargacc')
      SFCVAR(tbargsacc, 'tbargsacc')
      SFCVAR(tcanoaccgat, 'tcanoaccgat')
      SFCVAR(tcansacc, 'tcansacc')
      SFCVAR(tcoldm, 'tcoldm')
      SFCVAR(tcurm, 'tcurm')
      SFCVAR(thicecacc, 'thicecacc')
      SFCVAR(thliqcacc, 'thliqcacc')
      SFCVAR(thliqgacc, 'thliqgacc')
      SFCVAR(tmonthb, 'tmonthb')
      SFCVAR(todfrac, 'todfrac')
      SFCVAR(twarmm, 'twarmm')
      SFCVAR(tymaxlai, 'tymaxlai')
      SFCVAR(uvaccgat, 'uvaccgat')
      SFCVAR(veghght, 'veghght')
      SFCVAR(vgbiomas, 'vgbiomas')
      SFCVAR(vvaccgat, 'vvaccgat')
      SFCVAR(wdmindex, 'wdmindex')
      SFCVAR(zolnc, 'zolnc')
      SFCVAR(afrleaf, 'afrleaf')
      SFCVAR(afrroot, 'afrroot')
      SFCVAR(afrstem, 'afrstem')
      SFCVAR(autores, 'autores')
      SFCVAR(autresveg, 'autresveg')
      SFCVAR(burnarea, 'burnarea')
      SFCVAR(colrate, 'colrate')
      SFCVAR(dstcemls, 'dstcemls')
      SFCVAR(dstcemls3, 'dstcemls3')
      SFCVAR(gpp, 'gpp')
      SFCVAR(gppveg, 'gppveg')
      SFCVAR(grclarea, 'grclarea')
      SFCVAR(hetresveg, 'hetresveg')
      SFCVAR(hetrores, 'hetrores')
      SFCVAR(humiftrs, 'humiftrs')
      SFCVAR(leaflitr, 'leaflitr')
      SFCVAR(litrfall, 'litrfall')
      SFCVAR(litres, 'litres')
      SFCVAR(litresveg, 'litresveg')
      SFCVAR(lucemcom, 'lucemcom')
      SFCVAR(lucltrin, 'lucltrin')
      SFCVAR(lucsocin, 'lucsocin')
      SFCVAR(ltstatus, 'ltstatus')
      SFCVAR(mortrate, 'mortrate')
      SFCVAR(nbp, 'nbp')
      SFCVAR(nbpveg, 'nbpveg')
      SFCVAR(nep, 'nep')
      SFCVAR(nepveg, 'nepveg')
      SFCVAR(npp, 'npp')
      SFCVAR(probfire, 'probfire')
      SFCVAR(rg, 'rg')
      SFCVAR(rgveg, 'rgveg')
      SFCVAR(rm, 'rm')
      SFCVAR(rml, 'rml')
      SFCVAR(rmlvegacc, 'rmlvegacc')
      SFCVAR(rmr, 'rmr')
      SFCVAR(rmrveg, 'rmrveg')
      SFCVAR(rms, 'rms')
      SFCVAR(rmsveg, 'rmsveg')
      SFCVAR(roottemp, 'roottemp')
      SFCVAR(socres, 'socres')
      SFCVAR(socresveg, 'socresveg')
      SFCVAR(soilresp, 'soilresp')
      SFCVAR(tltrleaf, 'tltrleaf')
      SFCVAR(tltrroot, 'tltrroot')
      SFCVAR(tltrstem, 'tltrstem')
      SFCVAR(vgbiomas_veg, 'vgbiomas_veg')
      SFCVAR(wtstatus, 'wtstatus')
      ! End of CTEM fields

   end type SFCVARLIST_T

   integer, parameter :: INDX_SOIL    =  1
   integer, parameter :: INDX_GLACIER =  2
   integer, parameter :: INDX_WATER   =  3
   integer, parameter :: INDX_ICE     =  4
   integer, parameter :: INDX_AGREGE  =  5
   integer, parameter :: INDX_URB     =  6
   integer, parameter :: INDX_LAKE    =  7
   integer, parameter :: INDX_MAX     =  7

   type(SFCVARLIST_T), target  :: vd
   type(SFCVAR_T), allocatable :: vl(:)
   type(sfcptr), allocatable :: busptr(:)
   integer, allocatable :: statut(:,:)

   integer :: surfesptot = 0
   integer :: nvarsurf = 0  !# Number of surface bus var
   integer :: nsurf    = 0  !# Number of surface "types"
   integer :: tsrad_i=0, z0_i=0, z0t_i=0 !#TODO: remove, replace by vd%tsrad...

   integer :: drain=0
   integer :: drainaf=0
   integer :: insmavg=0
   integer :: isoil=0
   integer :: leg=0
   integer :: legaf=0
   integer :: ler=0
   integer :: leraf=0
   integer :: les=0
   integer :: lesaf=0
   integer :: letr=0
   integer :: letraf=0
   integer :: lev=0
   integer :: levaf=0
   integer :: overfl=0
   integer :: overflaf=0
   integer :: rootdp=0
   integer :: wflux=0
   integer :: wfluxaf=0
   integer :: wsoil=0

contains


   function sfcbus_init() result(F_istat)
      use phy_typedef, only: phymeta
      use phygetmetaplus_mod, only: phymetaplus, phygetmetaplus
      use sfc_options, only: schmurb, schmlake
      implicit none
      integer :: F_istat

#include <msg.h>
#include <rmnlib_basics.hf>
#include <clib_interface_mu.hf>

      integer :: i, istat, mulmax, idxmax
      type(SFCVAR_T) :: vl0(1)
      type(phymeta) :: mymeta
      type(phymetaplus) :: mymetaplus

      F_istat = RMN_ERR

      if (nsurf == 0) then
         idxmax = max(INDX_SOIL, INDX_GLACIER, INDX_WATER, INDX_ICE, INDX_AGREGE)
         if (schmurb  /= 'NIL') idxmax = max(idxmax, INDX_URB )
         if (schmlake /= 'NIL') idxmax = max(idxmax, INDX_LAKE)
         nsurf = idxmax - 1
      endif
 
      nvarsurf = size(transfer(vd, vl0))
      allocate(vl(nvarsurf))
      allocate(busptr(nvarsurf))
      vl = transfer(vd, vl)
      mulmax = 0
      do i = 1,nvarsurf
         vl(i)%i = i
         istat = clib_toupper(vl(i)%n)
         nullify(busptr(i)%ptr)
         istat = phygetmetaplus(mymetaplus, vl(i)%n, F_npath='V', &
              F_bpath='DPVE', F_quiet=.true., F_shortmatch=.false.)
         if (istat >= 0) then
            mymeta = mymetaplus%meta
            busptr(i)%ptr => mymetaplus%vptr
            vl(i)%doagg_L = (mymeta%bus(1:1) /= 'E')
            vl(i)%mul = mymeta%fmul
            vl(i)%niveaux = mymeta%nk
            vl(i)%mosaik = mymeta%mosaic + 1
            mulmax = max(mulmax, vl(i)%mul)
         endif

         select case(vl(i)%n)
         case('TSRAD')
            tsrad_i = i
         case('Z0')
            z0_i = i
         case('Z0T')
            z0t_i = i
         case('DRAIN')
            drain = mymetaplus%index
         case('DRAINAF')
            drainaf = mymetaplus%index
         case('INSMAVG')
            insmavg = mymetaplus%index
         case('ISOIL')
            isoil = mymetaplus%index
         case('LEG')
            leg = mymetaplus%index
         case('LEGAF')
            legaf = mymetaplus%index
         case('LER')
            ler = mymetaplus%index
         case('LERAF')
            leraf = mymetaplus%index
         case('LES')
            les = mymetaplus%index
         case('LESAF')
            lesaf = mymetaplus%index
         case('LETR')
            letr = mymetaplus%index
         case('LETRAF')
            letraf = mymetaplus%index
         case('LEV')
            lev = mymetaplus%index
         case('LEVAF')
            levaf = mymetaplus%index
         case('OVERFL')
            overfl = mymetaplus%index
         case('OVERFLAF')
            overflaf = mymetaplus%index
         case('ROOTDP')
            rootdp = mymetaplus%index
         case('WFLUX')
            wflux = mymetaplus%index
         case('WFLUXAF')
            wfluxaf = mymetaplus%index
         case('WSOIL')
            wsoil = mymetaplus%index
         end select

      enddo
      vd = transfer(vl, vd)
 
      allocate(statut(nvarsurf, mulmax))
      statut = 0

      F_istat = RMN_OK
      return
   end function sfcbus_init

end module sfcbus_mod
