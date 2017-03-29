    
pro makejeans_p12,mbh
	    file='2mass_galfit/PGC012557_aK_asky_981029n1280021_galfit.fits'
	    sersics=read_galfit_sersics(file)
	    Msun=3.27 ; Solar luminsoity in Ks band
	    a_b=0.05 ; K from schlafly
;		sersics.q = sersics.q>0.9
	    mge=sersics2mge(sersics,msun,a_b=a_b)
	    distance=67.0 ; mpc
	    inc=89.0
	    normpsf=[1.0]
	    sigmapsf=[0.16]
	    ;normpsf=[0.6,0.4]
	    ;sigmapsf=[0.07,0.4]
		;stop
		betaz=0.0
        ; Gauss+Moffat: 31% of the light in gaussian of 0.135” FWHM, 69% of the light in Moffat profile of FWHM=1.80", index fixed=4.765.  See attached figure for fit.
		normpsf=[0.31,1-0.31]
		sigmapsf=[0.135,1.8]/2.35

        readcol,'nifs_kinematics/pgc12557_wallace_init_sn40.dat',BinNum, X, Y, Npix,SN_expected,SN_measured,chi2,vbin,sbin,h3,h4
;        a=mrdfits('nifs_kinematics/NGC0722_binning_stellar_kinematics_wingeGNIRS_LSF_noemission_2moments_error_tSN55_deg4_1arcsec.fits',1)

        rotate_points, x, y, -sersics[0].pa+15, xbar, ybar		

		er_vbin=vbin*0+10
		er_sbin=sbin*0+10
		
	    pixsize=0.05
;		pixang=-57
		vbin=vbin-median(vbin)
		Vrms = sqrt(vbin^2 + sbin^2) 
		ERMS = sqrt((er_vbin*vbin)^2 + (er_sbin*sbin)^2)/VRMS
		plot_velfield,xbar,ybar,vbin,range=[-100,100]
		plot_velfield,xbar,ybar,vrms
		
		r=sqrt(xbar^2+ybar^2)
		; mark the central bin as off. The jeans model get the dispersion wrong here.
		;c=min(r,s)
		;erms[s] = sqrt(erms[s] + 160.0^2)
;		xbar[s] += 0.01
	    mbh1=1e6   
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh1, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf,ml=ml1
        
			!p.multi=[0,1,1]
			
	    plot,r,vrms,psym=3,/yno,xtitle='radius (arcsec)',ytitle='VRMS (km/s)',title='PGC12557',xrange=[0,1.3];,yrange=[100,160]
	    oploterror,r,vrms,psym=3,/nohat,erms
		oplot,r,rmsmodel,psym=1,color=127
	    print,mbh1*ml1,ml1,total(vrms-rmsmodel)^2/total(erms^2)

;		mbh= 1e5
		chi = 1e300
	    for lmbh=7.0,10,0.1 do begin
		mbh=10d^float(lmbh    )
	    jam_axisymmetric_rms, $
	        mge[*,0], mge[*,1], mge[*,2], mge[*,0], mge[*,1], mge[*,2], $
	        inc, mbh, distance, xbar, ybar, rmsModel, $
	        BETA=betaz, RMS=vrms, erms=erms,SIGMAPSF=psf,flux=flux, PIXSIZE=pixsize,pixang=pixang,normpsf=normpsf,ml=mlo
            ; stop iterating once we hit a detoration of 1% is chi^2
			print,chi,total(vrms-rmsmodel)^2/total(erms^2),alog10(mbh*mlo)
			if chi*1.01 lt total(vrms-rmsmodel)^2/total(erms^2) then break
			chi = min([total(vrms-rmsmodel)^2/total(erms^2),chi])
			ml2=temporary(mlo)
			rmsmodel2=rmsmodel
			mbh2=mbh
	endfor

        
	    oplot,r,rmsmodel,psym=4,color=200

		al_legend,['bh='+strn(alog10(mbh1*ml1),format='(G5.2)')+' log(Msun)','bh='+strn(alog10(mbh2*ml2),format='(G5.2)')+' log(Msun)'],/right,psym=[1,4]
    WRITE_PNG, 'PGC012557_jeans.png', TVRD(/TRUE)
	    stop
	end
    
    