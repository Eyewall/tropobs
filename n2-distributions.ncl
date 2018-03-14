;---------------------------------------------
; n2-distributions.ncl
; Patrick Duran
; 02 February 2017
; 
; Pulls data from a netCDF file containing
; dropsonde or rawinsonde data and computes
; the squared Brunt-Vaisala frequency.
; Plots distributions of N^2 within a
; selected vertical layer.
; Computes a sharpness parameter and also
; plots distributions of that.
;
; 06 March 2017: Added composite averaging with respect to the
;                cold-point tropopause height.
;---------------------------------------------

infile = "/nfs/jm13/pduran/rawinsondes/sparc/rawinsondes_0-1000km.nc"
;infile = "/home/patrick/research/data/rawinsonde/rawinsondes_0-1000km.nc"
cat    = "hur"  ;Storm category for which to compute averages
                ;Options: "td", "ts", "td+ts", "hur","cat1+2","maj","all"
avgtype= 1      ;0=Plot simple average; 1=Plot composite average WRT CPT height
outdir = "/nfs/jm13/pduran/web/trop/"+cat+"/"
;outdir = "/home/patrick/research/plots/trop/"+cat+"/"
minrad = 0      ;Minimum radius outside of which to pull sonde data
maxrad = 1000    ;Maximum radius within which to pull sonde data
dr     = 100    ;Width of radial bands
mincpz = 14000  ;Minimum cold-point tropopause height allowed (m)
g = 9.80665
rads  = ispan(minrad,maxrad,dr)
nrads = (maxrad-minrad)/dr

system("mkdir "+outdir)

;Read input file
print( "Reading input file "+infile+"   "+systemfunc("date"))
ncfile = addfile(infile,"r")
droptimein = ncfile->droptime
droplatin  = ncfile->droplat
droplonin  = ncfile->droplon
heightin   = ncfile->height
cptin      = ncfile->cpt
cpzin      = ncfile->cpz
tophghtin  = ncfile->maxhght
presin     = ncfile->p
tcin       = ncfile->t
sondelatin = ncfile->sondelat
sondelonin = ncfile->sondelon
radiusin   = ncfile->radius
azimuthin  = ncfile->azimuth
stormlatin = ncfile->stormlat
stormlonin = ncfile->stormlon
vmaxin     = ncfile->vmax
pminin     = ncfile->pmin

;Total number of sondes in the entire dataset
totnumsondes = dimsizes(droptimein)

;Compute squared Brunt-Vaisala frequency
tk = tcin+273.15
th = pot_temp(presin,tk,-1,False)
dthdz = center_finite_diff_n(th,heightin,False,0,1)  ;Units K/m
n2in = (g/th)*dthdz*10000 ;Units 10^-4 s^-2

;If the cold-point height is smaller than the user-defined
;minimum value, set it and the cold point temperature to missing.
cpzin = where(cpzin.lt.mincpz,cpzin@_FillValue,cpzin)
cptin = where(cpzin.lt.mincpz,cptin@_FillValue,cptin)

;If the cold point is at the very top of the sounding, set it
;to missing.
cpzin = where(cpzin.eq.tophghtin,cpzin@_FillValue,cpzin)
cptin = where(cpzin.eq.tophghtin,cptin@_FillValue,cptin)

;Vertical grid spacing and midpoint of vertical coordinate
dz = heightin(1)-heightin(0)
midz = heightin(dimsizes(heightin)/2)
nhghts  = dimsizes(heightin)

;Initialize shifted arrays for composite averaging
if(avgtype.eq.1)then
   tcshift = new(nhghts*2-1,"float")
   tcshift!0 = "comphght"
   tcshift&comphght = ispan(max(heightin)*-1,max(heightin),dz)
   n2shift = tcshift
end if

;Initialize average arrays
pres_avg= new((/nrads,nhghts/),"float")
tc_avg  = new((/nrads,nhghts/),"float")
n2_avg  = new((/nrads,nhghts/),"float")
cpz_cat = new((/nrads,totnumsondes/),"float")
cpt_cat = new((/nrads,totnumsondes/),"float")
nsondes = new(nrads,"integer")
cpz_avg = new(nrads,"float")
cpt_avg = new(nrads,"float")

;Get data only within the selected radial band
do i=0,nrads-1
   minbinrad = rads(i)
   maxbinrad = rads(i+1)
   radinds  := ind(radiusin.ge.minbinrad.and.radiusin.lt.maxbinrad)

   droptime := droptimein(radinds) 
   droplat  := droplatin (radinds) 
   droplon  := droplonin (radinds) 
   height   := heightin 
   cpt      := cptin     (radinds) 
   cpz      := cpzin     (radinds) 
   tophght  := tophghtin (radinds) 
   pres     := presin    (radinds,:) 
   tc       := tcin      (radinds,:) 
   sondelat := sondelatin(radinds,:) 
   sondelon := sondelonin(radinds,:) 
   radius   := radiusin  (radinds) 
   azimuth  := azimuthin (radinds) 
   stormlat := stormlatin(radinds) 
   stormlon := stormlonin(radinds) 
   vmax     := vmaxin    (radinds) 
   pmin     := pminin    (radinds) 
   n2       := n2in      (radinds,:)

   n2!0 = "radius"
   n2!1 = "height"
   n2&height = height

   ;Get indices for each storm category
   if(cat.eq."td")then
      catinds := ind(vmax.le.33)
   else if(cat.eq."ts")then
      catinds := ind(vmax.ge.34.and.vmax.lt.64)
   else if(cat.eq."td+ts")then
      catinds := ind(vmax.lt.64)
   else if(cat.eq."hur")then
      catinds := ind(vmax.ge.64)
   else if(cat.eq."cat1+2")then
      catinds := ind(vmax.ge.64.and.vmax.le.95)
   else if(cat.eq."maj")then
      catinds := ind(vmax.gt.95)
   else if(cat.eq."all")then
      catinds := ind(vmax.ge.0)
   end if
   end if
   end if
   end if
   end if
   end if
   end if

   ;Get cold-point tropopause temperatures and heights for each sonde
   ;in this radial bin and category.
   cpz_cat(i,:dimsizes(catinds)-1) = cpz(catinds)
   cpt_cat(i,:dimsizes(catinds)-1) = cpt(catinds)
   ;Get the number of sondes with data up to the cold point or higher
   nsondes(i)= dimsizes(ind(.not.ismissing(cpz_cat(i,:))))
   ;Compute some mean fields for this radial band
   ;and category
   if(avgtype.eq.0)then ;Compute simple average at every grid point
      pres_avg(i,:) = dim_avg_n_Wrap(pres(catinds,:),0)
      tc_avg  (i,:) = dim_avg_n_Wrap(tc  (catinds,:),0)
      n2_avg  (i,:) = dim_avg_n_Wrap(n2  (catinds,:),0)
   else if(avgtype.eq.1) ;Compute composite averge with respect to cold-point height
      do j=0,dimsizes(catinds)-1
         catind = catinds(j)
         cpzcur = cpz_cat(i,j)      ;Current cold point height
         if(.not.ismissing(cpzcur))
            print( cpzcur)
            maxhghtshift = max(height)-cpzcur
            minhghtshift = -cpzcur+dz
            tccur = tc(catind,:)    ;T profile for this sonde
            n2cur = n2(catind,:)    ;N^2 profile for this sonde
            tcshift({0:maxhghtshift}) = (/tccur({cpzcur:})/)
            tcshift({minhghtshift:0}) = (/tccur({:cpzcur-dz})/)
            n2shift({0:maxhghtshift}) = (/n2cur({cpzcur:})/)
            n2shift({minhghtshift:0}) = (/n2cur({:cpzcur-dz})/)
            print(n2shift)
            exit
         end if
      end do
      exit
   end if
   end if

   ;Determine cold point tropopause temperature and height using mean fields
   cpz_avg(i) = height(minind(tc_avg(i,:)))/1000. ;Convert to km
   cpt_avg(i) = min(tc_avg(i,:))
end do

;Assign coordinate variables
midpoints = new(dimsizes(rads)-1,"integer")
do i=0,nrads-1
   midpoints(i) = (rads(i+1)+rads(i))/2
end do
nsondes!0 = "radius"
nsondes&radius = midpoints
copy_VarCoords(nsondes,cpz_avg)
copy_VarCoords(nsondes,cpt_avg)
n2_avg!0 = "radius"
n2_avg!1 = "height"
n2_avg&radius = midpoints
n2_avg&height = heightkm

;Compute probability distributions
cpzopt = True
cpzopt@bin_min = 14.
cpzopt@bin_max = 19.
cpzopt@bin_nice= True
pdf_cpz = pdfx(cpz_cat,21,cpzopt)
cptopt = cpzopt
cptopt@bin_min = -90
cptopt@bin_max = -50
pdf_cpt = pdfx(cpt_cat,61,cptopt)

;-------------------;
;Cross-section plots;
;-------------------;
out_suffix = minrad+"-"+maxrad+"km"
;STATIC STABILITY AND COLD-POINT TROPOPAUSE HEIGHT
;Transpose array so that height is leftmost dimension
n2_avg_plot = transpose(n2_avg)
wks_stab = gsn_open_wks("png",outdir+"stab+cp_trop_hght_"+out_suffix)
gsn_define_colormap(wks_stab,"WhiteGreen")
stabres = True
stabres@cnFillOn = True
stabres@cnLevelSelectionMode = "ManualLevels"
stabres@cnMinLevelValF = 0
stabres@cnMaxLevelValF = 7.0
stabres@cnLevelSpacingF = 1.0
stabres@cnLineLabelsOn = False
stabres@cnLinesOn = False
stabres@gsnLeftStringFontHeightF = 0.014
stabres@gsnRightStringFontHeightF = stabres@gsnLeftStringFontHeightF
stabres@gsnLeftString = "N~S~2~N~ (10~S~-4~N~ s~S~-2~N~) and cold point height"
stabres@gsnRightString = "Category: "+str_upper(cat)+"        "+minrad+"-"+maxrad+" km"
stabres@tiXAxisString = "Radius (km)"
stabres@tiYAxisString = "Height (km)"
stabres@lbLabelBarOn = True
stabres@lbAutoManage = False
stabres@lbOrientation = "Vertical"
stabres@cnInfoLabelOn = False
stabres@trYMinF = 9
stabres@trYMaxF = 25
stabres@tmYROn = False
stabres@tmXTOn = False
stabres@gsnMaximize = True
stabres@gsnDraw = False
stabres@gsnFrame = False
stab_img = gsn_csm_contour(wks_stab,n2_avg_plot,stabres)
tropres = True
tropres@gsLineColor = "yellow"
tropres@gsLineThicknessF = 18
trophghts = gsn_add_polyline(wks_stab,stab_img,cpz_avg&radius,cpz_avg,tropres)
draw(stab_img)
frame(wks_stab)

;-------------;
;Distributions;
;-------------;
;COLD POINT TROPOPAUSE HEIGHT
pdf_cpz_wks = gsn_open_wks("png",outdir+"pdf_cpz_"+out_suffix)
pdf_cpt_wks = gsn_open_wks("png",outdir+"pdf_cpt_"+out_suffix)
pdfres = True
;pdfres@tmYROn = False
;pdfres@tmXTOn = False
pdfres@gsnMaximize = True
pdfres@tiXAxisOn = False
pdfres@xyLineThicknessF = 6.
pdfres@tiMainFontHeightF = 0.02
pdfres@tiMainString = "PDF of Cold-Point Tropopause Height for "+str_upper(cat)
pdfres@tiYAxisFontHeightF = 0.02
pdfres@tiXAxisFontHeightF = pdfres@tiYAxisFontHeightF
pdfres@tiYAxisString = "Percent"
pdfres@tiXAxisString = "Cold-Point Tropopause Height (km)"
pdfres@trXMinF = cpzopt@bin_min;pdf_cpz@bin_bound_min
pdfres@trXMaxF = cpzopt@bin_max;pdf_cpz@bin_bound_max
pdfres@trYMaxF = 30
pdf_cpz_plot = gsn_csm_xy(pdf_cpz_wks,pdf_cpz@bin_center,pdf_cpz,pdfres)
pdfres@tiMainString = "PDF of Cold-Point Tropopause Temperature for "+str_upper(cat)
pdfres@tiXAxisString = "Cold-Point Tropopause Temperature (degrees C)"
pdfres@trXMinF = cptopt@bin_min;pdf_cpz@bin_bound_min
pdfres@trXMaxF = cptopt@bin_max;pdf_cpz@bin_bound_max
pdfres@trYMaxF = 11
pdf_cpt_plot = gsn_csm_xy(pdf_cpt_wks,pdf_cpt@bin_center,pdf_cpt,pdfres)

;print( nsondes)
;printVarSummary( n2_avg_plot)
;print( cpz_avg)
;print( cpt_avg)
;print( sum(nsondes))