    
pro makejeans_n1222,mbh
	    wset,0
	    file='2mass_galfit/NGC1022_aK_asky_981010s1060044_galfit.fits'
	    sersics=read_galfit_sersics(file)
	    Msun=3.27 ; Solar luminsoity in Ks band
	    a_b=0.169 ; K from schlafly
;		sersics.q = sersics.q>0.9
	    mge=sersics2mge(sersics,msun,a_b=a_b)
	    distance=19.9000 ; mpc
	    inc=40.0
	    normpsf=[1.0]
	    sigmapsf=[0.16]
	    normpsf=[0.6,0.4]
	    sigmapsf=[0.07,0.4]
		;stop
		betaz=0.0
        ; Gauss+Moffat: 31% of the light in gaussian of 0.135” FWHM, 69% of the light in Moffat profile of FWHM=1.80", index fixed=4.765.  See attached figure for fit.
		;normpsf=[0.31,1-0.31]
		;sigmapsf=[0.135,1.8]/2.35


  readcol,'nifs_kinematics/NGC1022_binning_stellar_kinematics_wingeGNIRS_LSF_noemission_secondMOM_errorSTD_tSN55_cut12_deg4.txt',$
      xBar1, yBar1, srn, vbin, sbin, er_vbin,er_sbin	
      rotate_points,xbar1,ybar1, 90-sersics[0].pa, xbar, ybar		
	  
;	    xbar+= 0.01
	    pixsize=0.05
;		pixang=-57
		vbin=vbin-median(vbin)
		Vrms = sqrt(vbin^2 + sbin^2) 
		ERMS = sqrt((er_vbin*vbin)^2 + (er_sbin*sbin)^2)/VRMS

		; mark the central bin as off. The jeans model get the dispersion wrong here.
		r=sqrt(xbar^2+ybar^2)
		c=min(r,s)
		erms[s] = sqrt(erms[s] + 160.0^2)

		plot_velfield,xbar,ybar,vbin,range=[-100,100]
		plot_velfield,xbar,ybar,vrms
		
	    mbh1=1e5   
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh1, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf,ml=ml1
        
			!p.multi=[0,1,1]
			r=sqrt(xbar^2+ybar^2)
	    plot,r,vrms,psym=3,/yno,xtitle='radius (arcsec)',ytitle='VRMS (km/s)',title='NGC1022',xrange=[0,1.5]
	    oploterror,r,vrms,psym=3,/nohat,erms
	    oplot,r,rmsmodel,psym=1,color=127
	    print,mbh1*ml1,ml1,total(vrms-rmsmodel)^2/total(erms^2)

		chi = 1e300
	    for lmbh=8.0,12,0.1 do begin
		mbh=10d^float(lmbh    )
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf,ml=mlo
            ; stop iterating once we hit a detoration of 1% is chi^2
			print,chi,total(vrms-rmsmodel)^2/total(erms^2),alog10(mbh*mlo),mbh
			if chi*1.01 lt total(vrms-rmsmodel)^2/total(erms^2) then break
			chi = min([total(vrms-rmsmodel)^2/total(erms^2),chi])
			ml2=temporary(mlo)
			rmsmodel2=rmsmodel
			mbh2=mbh
	endfor
    oplot,r,rmsmodel,psym=4,color=200
	
;	    oplot,r,rmsmodel,psym=4,color=200

;		mbh= 1.5848925e+08   
;;	    for lmbh=7.0,9,0.1 do begin
;;		mbh=10d^float(lmbh    )
;		print,mbh
;	    jam_axisymmetric_rms, $
;	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
;	        inc, mbh, distance, xbar, ybar, rmsModel, $
;	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf,ml=ml
;	    oplot,r,rmsmodel,psym=4,color=200
;;	endfor

    
	al_legend,['bh='+strn(alog10(1e5*ml1),format='(G5.2)')+' log(Msun)','bh='+strn(alog10(mbh2*ml2),format='(G5.2)')+' log(Msun)'],/right,psym=[1,4]
	WRITE_PNG, 'NGC1022_jeans.png', TVRD(/TRUE)
	    stop
	end
    
    