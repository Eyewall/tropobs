load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;---------------------------------------------
; tvsz.ncl
; Patrick Duran
; 10 April 2018
; 
; Pulls data from a netCDF file containing
; dropsonde or rawinsonde data and plots
; vertical profiles of temperature, the
; vertical temperature gradient, and the
; second derivative of temperature WRT height.
; Averages over all radial bins to get a single
; vertical profile for each quantity.
;
;---------------------------------------------


machine = "personal"  ;Set to "personal" or "albany"
cattype = "categories"        ;Set to "categories" or "tb"
if(cattype.eq."tb")then
   cats    = (/"coldest","cold","warm"/)
else
   cats    = (/"hur","td+ts"/)
end if
infiles = new(dimsizes(cats),"string")
do i=0,dimsizes(cats)-1
   if(machine.eq."personal")then
      infiles(i) = "/home/patrick/research/data/rawinsonde/rawinsondes_0-1000km_cpzcomposite_"+cats(i)+".nc" 
      outdir = "/home/patrick/research/plots/"
   else
      infiles(i) = "/jm13/pduran/rawinsondes/sparc/rawinsondes_0-1000km_cpzcomposite_"+cats(i)+".nc" 
      outdir = "/jm13/pduran/web/"
   end if
end do

;Read data
ncfiles = addfiles(infiles,"r")
ListSetType(ncfiles,"join")
t = ncfiles[:]->tc_avg
nsondes = ncfiles[:]->nsondes

;Compute vertical derivatives
dtdz   = center_finite_diff_n(t,t&comphght,False,0,2)
dt2dz2 = center_finite_diff_n(dtdz,t&comphght,False,0,2)
dtdz!0 = "cats"
dtdz!1 = "radius"
dtdz!2 = "comphght"
dtdz&cats     = cats
dtdz&radius   = t&radius
dtdz&comphght = t&comphght
copy_VarCoords(dtdz,dt2dz2)

;Average over radius dimension and sum the sondes
tavg      = dim_avg_n_Wrap(t,1)
dtdzavg   = dim_avg_n_Wrap(dtdz,1)
dt2dz2avg = dim_avg_n_Wrap(dt2dz2,1)
nsondestot= dim_sum(nsondes)

;Get strings we will use for legend labels
tbs = (/" T~B~b~N~ ~F34~#~F21~ -50~S~o~N~C",\
        " -50~S~o~N~C < T~B~b~N~ ~F34~#~F21~ -30~S~o~N~C",\
        " T~B~b~N~ > -30~S~o~N~C"/)
if(cattype.eq."tb")then
   catswithnumsondes = tbs+" (N="+tostring(nsondestot)+")"
else
   catswithnumsondes = str_upper(cats)+" (N="+tostring(nsondestot)+")"
end if

;Plot the vertical profiles
;Temperature
twks = gsn_open_wks("png",outdir+"tavg_"+str_concat(cats))
res= True
res@tiYAxisString = "Height above or below tropopause (km)"
res@tiXAxisString = "Temperature (~S~o~N~C)"
res@trXMinF = -76
res@trXMaxF = -58
res@trYMinF = -3
res@trYMaxF = 3
if(dimsizes(cats).eq.3)then
   res@xyLineColors = (/"magenta","red","blue"/)
   res@xyLineThicknesses = (/6,6,6/)*1.
   res@xyDashPatterns = (/0,0,0/)
else
   res@xyLineColors = (/"red","blue"/)
   res@xyLineThicknesses = (/6,6/)*1.
   res@xyDashPatterns = (/0,0/)
end if   
res@pmLegendDisplayMode = "Always"
res@pmLegendWidthF = 0.15
res@pmLegendHeightF = 0.15
res@pmLegendParallelPosF = 0.66
res@pmLegendOrthogonalPosF = -0.80
res@lgAutoManage = False
res@lgLabelFontHeightF = 0.02
res@lgPerimOn = False
res@lgBoxMinorExtentF = 0.3
res@xyExplicitLegendLabels = catswithnumsondes
res@gsnMaximize = True
tplot = gsn_csm_xy(twks,tavg,tavg&comphght,res)

;Vertical temperature gradient
dtdzwks = gsn_open_wks("png",outdir+"dtdzavg_"+str_concat(cats))
res@trXMinF = -10
res@trXMaxF = 10
res@tiXAxisString = "Vertical temperature gradient (~S~o~N~C km~S~-1~N~)"
res@pmLegendDisplayMode = "Never"
dtdzplot = gsn_csm_xy(dtdzwks,dtdzavg,dtdzavg&comphght,res)

;Second derivative of temperature with respect to height
dt2dz2wks = gsn_open_wks("png",outdir+"dt2dz2avg_"+str_concat(cats))
res@trXMinF = -20
res@trXMaxF = 60
res@tiXAxisString = "2nd derivative of T WRT height (~S~o~N~C km~S~-2~N~)"
dt2dz2plot = gsn_csm_xy(dt2dz2wks,dt2dz2avg,dt2dz2avg&comphght,res)
