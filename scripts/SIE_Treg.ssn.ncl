begin
dirSrfc = "/Users/irudeva/work/DATA/ERAint/Mslp_1deg/"

diro = "../output/SIE/" 

syr = 1980
eyr = 2017
nyr = eyr-syr+1
yrs = ispan(syr,eyr,1)

ssn = "JAS"
if (ssn.eq."JFM") then
   mStrt = 1 ; -month start
   mLast = 3; -month last
   maxndays = 90
   nmon = 3
   mon = (/1,2,3/)
   mday = (/31,28,31/)
else if  (ssn.eq."JAS") then
   mStrt = 7 ; -month start
   mLast = 9 ; -month last
   maxndays = 92
   nmon = 3
   mon = (/7,8,9/)
   mday = (/31,31,30/)
end if
end if

; regions

reg = (/"C_Europe ","E_Europe ","W_Siberia","E_Siberia","W_China  ","E_China  ","W_America","E_America","Kara_BarS","ESib_Chuk" /)
nreg = dimsizes(reg)

lonreg1 = (/ 5, 30, 70, 100,  80, 100, 235, 260, 30, 160/)
lonreg2 = (/25, 50, 90, 120, 100, 120, 255, 280, 70, 200/)
latreg1 = (/40, 45, 50,  55,  25,  22,  45,  35, 70,  70/)
latreg2 = (/55, 60, 65,  70,  45,  42,  60,  50, 80,  80/)

Tssn = new((/nyr,nreg+1/),float,-999)   ; spacially averaged temperature for each region
Tssn@_FillValue = -999

ccr      = new(nreg,float)
ccrdt    = ccr
ccrsig   = new(nreg,integer)
ccrdtsig = new(nreg,integer)

psig = 0.1  
; *********************************************
; read monthly SIE
; *********************************************
ncol = 4
do im = 0,dimsizes(mon)-1

        if (mon(im).lt.10)then
            month = "0"+mon(im)
        else
            month = ""+mon(im)
        end if
        print(month)
        data = readAsciiTable("../data/SIE/N_"+month+"_extent_v3.0.csv", ncol, "float", (/1,0/))

        ;create an empty SIE array
        if(im.eq.0)then
            monSIE =new((/dimsizes(mon),dimsizes(data(:,0))/),float)
            monSIE@_FillValue = -9999
            yearSIE = data(:,0)
        end if

        monSIE(im,:) = data(:,2)
end do

; select years
iy1 = ind(yearSIE.eq.syr)
iy2 = ind(yearSIE.eq.eyr)

; seasonal mean
SIEssn = dim_avg_n(monSIE(:,iy1:iy2), 0)


; *********************************************
; Temperature
; *********************************************
do yr = syr,eyr

  fname = systemfunc("ls "+dirSrfc+"erain.mslp_air."+ssn+"mean."+yr+".nc")
  print ("reading file "+ fname)
  clfile = addfile(fname,"r")

  if (getfilevartypes(clfile,"t2m").eq."short")then
   tvar = short2flt(clfile->t2m)
  else
   tvar = clfile->t2m
  end if
  tvar = tvar - 273.15
  tvar@units = "degC" 

  lon  = clfile->lon
  lat  = clfile->lat

  do ireg = 0, nreg-1
    treg = tvar(0,{latreg1(ireg):latreg2(ireg)},{lonreg1(ireg):lonreg2(ireg)})
    Tssn(yr-syr,ireg) =avg(treg)

    delete(treg)
  end do

end do


; *********************************************
; trends
; *********************************************

Tssn(:,nreg) = SIEssn ; adding SIE to Treg for trends
; detrending
dtTssn = dtrend_msg_n(yrs-syr,Tssn,False,True,0)
dtSIEssn  = dtrend_msg(yrs-syr,SIEssn,False,True)

Tslope   = dtTssn@slope
SIEslope = dtSIEssn@slope


;evaluating significance
tval = new(nreg+1,"float")  ; preallocate tval as a float array and
df   = new(nreg+1,"integer")  ; df as an integer array for use in regcoef 

rc = regcoef(yrs-syr,transpose(Tssn),tval,df)   ; regress z against a straight line to
                                                ; return the tval and degrees of freedom
df = equiv_sample_size(transpose(Tssn),0.05,0)  ; If your data may be significantly autocorrelated
                                    ; it is best to take that into account, and one can
                            ; do that by using equiv_sample_size. Note that 
                                    ; in this example df (output from regcoef) is 
                                    ; overwritten with the output from equiv_sample_size.
                                    ; If your data is not significantly autocorrelated one
                                    ; can skip using equiv_sample_size.

df = df-2          ; regcoef/equiv_sample_size return N, need N-2
beta_b = new(nreg+1,"float")    ; preallocate space for beta_b
beta_b = 0.5       ; set entire beta_b array to 0.5, the suggested value of beta_b 
                    ; according to betainc documentation
z_signif = (1.-betainc(df/(df+tval^2), df/2.0, beta_b))*100. ; significance of trends 
                                                            ; expressed from 0 to 100%

print(z_signif)


; *********************************************
; Correlations
; *********************************************
do ireg = 0,nreg-1
    ccr(ireg) = escorc(Tssn(:,ireg),SIEssn)
    ; test significance level    
    t    = ccr(ireg)*sqrt((nyr-2)/(1-ccr(ireg)^2))   
    p    = student_t(t, nyr-2)
    if (p.le.psig) then
        ccrsig(ireg) = 1
        ;print("r="+ccr(ireg)+" is significant at the 95% level")
    else
        ccrsig(ireg) = 0
        ;print("r="+ccr(ireg)+" is NOT significant at the 95% level")
    end if
end do

do ireg = 0,nreg-1
    ccrdt(ireg) = escorc(dtTssn(:,ireg),dtSIEssn)
    ; test significance level    
    t    = ccrdt(ireg)*sqrt((nyr-2)/(1-ccrdt(ireg)^2))   
    p    = student_t(t, nyr-2)
    if (p.le.psig) then
        ccrdtsig(ireg) = 1
        ;print("r="+ccr(ireg)+" is significant at the 95% level")
    else
        ccrdtsig(ireg) = 0
        ;print("r="+ccr(ireg)+" is NOT significant at the 95% level")
    end if
end do

; print(ccr)
; print(ccrdt)


; *********************************************
; write an ascii file
; *********************************************

hlist = [/"reg", "corr","corrsig","corr_dt","corr_dtsig","trend","trend_sig"/]

alist = [/reg, ccr, ccrsig, ccrdt, ccrdtsig, Tslope(:nreg-1), z_signif(:nreg-1)/]

fname = "../output/SIE/corr.ssnSIR_Treg."+ssn+".txt"
write_table(fname, "w", hlist, "%s%13s%8s%7s%7s%10s%s%s")
write_table(fname, "a", alist, "%s%7.2f%8d%7.2f%10d%7.3f%7.1f")


end