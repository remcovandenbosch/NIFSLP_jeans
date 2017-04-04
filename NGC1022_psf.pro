    
pro makejeans_n1022,mbh

	    file='2mass_galfit/NGC1022_aK_asky_981010s1060044_galfit.fits'
	    sersics=read_galfit_sersics(file)
	    Msun=3.27 ; Solar luminsoity in Ks band
	    a_b=0.169 ; K from schlafly
;		sersics.q = sersics.q>0.9
	    mge=sersics2mge(sersics,msun,a_b=a_b)
	    distance=19.9000 ; mpc
	    inc=40.0
	    normpsf1=[1.0]
	    sigmapsf1=[0.16]
	    ;normpsf=[0.6,0.4]
	    ;sigmapsf=[0.07,0.4]
		betaz=0.0
        ; Gauss+Moffat: 31% of the light in gaussian of 0.135” FWHM, 69% of the light in Moffat profile of FWHM=1.80", index fixed=4.765.  See attached figure for fit.
		normpsf2=[0.31,1-0.31]
		sigmapsf2=[0.135,1.8]/2.35


  readcol,'nifs_kinematics/NGC1022_binning_stellar_kinematics_wingeGNIRS_LSF_noemission_secondMOM_errorSTD_tSN55_cut12_deg4.txt',$
      xBar1, yBar1, srn, vbin, sbin, er_vbin,er_sbin	
      rotate_points,xbar1,ybar1, 90-sersics[0].pa, xbar, ybar		
	  
;	    xbar+= 0.01
	    pixsize=0.05
		pixang=-57.0
		vbin=vbin-median(vbin)
		Vrms = sqrt(vbin^2 + sbin^2) 
		ERMS = sqrt((er_vbin*vbin)^2 + (er_sbin*sbin)^2)/VRMS

		; mark the central bin as off. The jeans model get the dispersion wrong here.
		r=sqrt(xbar^2+ybar^2)
		c=min(r,s)
		erms[s] = sqrt(erms[s] + 160.0^2)

		plot_velfield,xbar,ybar,vbin,range=[-100,100]
		plot_velfield,xbar,ybar,vrms
		
;	    mbh1=1e5   
		chi = 1e300
	    for lmbh=8.0,12,0.1 do begin
		mbh=10d^float(lmbh    )
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf1,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf1,ml=mlo
            ; stop iterating once we hit a detoration of 1% is chi^2
			chinew=total((vrms-rmsmodel)^2/(erms^2))
			print,chi,chinew,alog10(mbh*mlo),mbh
			if chi*1.01 lt chinew then break
			chi = min([chinew,chi])
			ml1=temporary(mlo)
			rmsmodel1=rmsmodel
			mbh1=mbh
	endfor
			!p.multi=[0,1,1]
			r=sqrt(xbar^2+ybar^2)
	    plot,r,vrms,psym=3,/yno,xtitle='radius (arcsec)',ytitle='VRMS (km/s)',title='NGC1022',xrange=[0,0.99]
	    oploterror,r,vrms,psym=3,/nohat,erms
	    oplot,r,rmsmodel1,psym=1,color=127
	    print,mbh1*ml1,ml1,total(vrms-rmsmodel)^2/total(erms^2)
		print,'==========='


		chi = 1e300
	    for lmbh=8.0,12,0.1 do begin
		mbh=10d^float(lmbh    )
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf2,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf2,ml=mlo
            ; stop iterating once we hit a detoration of 1% is chi^2
			chinew=total((vrms-rmsmodel)^2/(erms^2))
			print,chi,chinew,alog10(mbh*mlo),mbh
			if chi*1.01 lt chinew then break
			chi = min([chinew,chi])
			ml2=temporary(mlo)
			rmsmodel2=rmsmodel
			mbh2=mbh
	endfor
    oplot,r,rmsmodel,psym=4,color=200
	
	al_legend,['Single Gauss PSF, bh='+strn(alog10(mbh1*ml1),format='(G5.2)')+' log(Msun)','Double Gauss PSF, bh='+strn(alog10(mbh2*ml2),format='(G5.2)')+' log(Msun)'],/right,psym=[1,4]
	WRITE_PNG, 'NGC1022_jeans_psf.png', TVRD(/TRUE)
	    stop
	end
    
    