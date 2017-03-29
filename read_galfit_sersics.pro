

function read_galfit_sersics,fitsfile,includeagn=includeagn,sky=sky

    model=mrdfits(fitsfile,2,hg)
    ;resid=mrdfits(name,3)

    ; find out how many sersics are in this model
    xc =float((strsplit(sxpar(hg,strn(2)+'_XC') ,'*+/-[]{}',/ex))[0])
    yc =float((strsplit(sxpar(hg,strn(2)+'_YC') ,'*+/-[]{}',/ex))[0])
    
    for i=1,1d4 do begin
	 if (xc-float((strsplit(sxpar(hg,strn(i+2)+'_XC') ,'*+/-[]{}',/ex))[0])) gt 0.1 then break
	 if (yc-float((strsplit(sxpar(hg,strn(i+2)+'_YC') ,'*+/-[]{}',/ex))[0])) gt 0.1 then break
	 ;if sxpar(hg,strn(i+2)+'_RE') eq 0 then break
	 if (strsplit(sxpar(hg,strn(i+2)+'_RE') ,'*+/-[]{}',/ex))[0] lt 1d-6 then break
    endfor
    nsersic=i

    if keyword_set(includeagn) then $
         sersics=replicate({mag:7.0,re:10.0,n:1.0,q:0.5,pa:0.0,x:256.0,y:512.0},nsersic+1) $
      else $
         sersics=replicate({mag:7.0,re:10.0,n:1.0,q:0.5,pa:0.0,x:256.0,y:512.0},nsersic)

	for i=0,nsersic-1 do begin
     sersics[i].mag =float((strsplit(sxpar(hg,strn(i+2)+'_MAG'),'*+/-[]{}',/ex))[0])
     sersics[i].re  =float((strsplit(sxpar(hg,strn(i+2)+'_RE') ,'*+/-[]{}',/ex))[0])
     sersics[i].n   =float((strsplit(sxpar(hg,strn(i+2)+'_N')  ,'*+/-[]{}',/ex))[0])
     sersics[i].q   =float((strsplit(sxpar(hg,strn(i+2)+'_AR') ,'*+/-[]{}',/ex))[0])
     sersics[i].pa  =float((strsplit(sxpar(hg,strn(i+2)+'_PA') ,'*+/[]{}',/ex))[0])
     sersics[i].x   =float((strsplit(sxpar(hg,strn(i+2)+'_XC') ,'*+/-[]{}',/ex))[0])
     sersics[i].y   =float((strsplit(sxpar(hg,strn(i+2)+'_YC') ,'*+/-[]{}',/ex))[0])
	  endfor 

     if keyword_set(includeagn) then begin
         agnmag=float((strsplit(sxpar(hg,strn(nsersic+2)+'_MAG') ,'*+/-[]{}',/ex))[0])
         sersics[nsersic].mag =agnmag[0]
         sersics[nsersic].re  =0.1
         sersics[nsersic].n   =0.5
         sersics[nsersic].q   = max(sersics[0:nsersic-1].q,c) ;0.98
         sersics[nsersic].pa  =sersics[c].pa
         sersics[nsersic].x   =sersics[c].x
         sersics[nsersic].y   =sersics[c].y
     endif

     ; find x and y offset relative to fitting box
     ;xoff=(float((strsplit(sxpar(hg,'FITSECT') ,'*+/-[]{},:',/ex))[0]))-1 
     ;yoff=(float((strsplit(sxpar(hg,'FITSECT') ,'*+/-[]{},:',/ex))[2])) -1

     sersics.x -= 1
     sersics.y -= 1

     sky=float((strsplit(sxpar(hg,'1_SKY') ,'*+/-[]{},:',/ex)) [0])
    
return, sersics
end

