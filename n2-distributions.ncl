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
; 05 April 2017: Added IBTrACS storm type to ensure that all storms
;                included here are tropical.
;---------------------------------------------


machine = "personal"  ;Set to "personal" or "albany"
intorir = 1          ;Set to 0 to get composite averages for storm categories
                     ;Set to 1 to get composite averages for IR brightness temps
cat     = "hur"      ;Storm category for which to compute averages
                     ;Options: "td", "ts", "td+ts", "hur","cat1+2","maj","all"
                     ;Used only if intorit=0
irscene = "coldest"  ;IR satellite scene for composite averaging
                     ;Options: "coldest", "cold", "warm"
                     ;Used only if intorir=1
coldcutoff = -50     ;IR brightness temperature bin cutoff temperature (deg C)
warmcutoff = -30     ;Used only if intorir=1
intype  = 0     ;0=Read rawinsonde input file; 1=Read composite average input file
type    = 1     ;0=Plot data; 1=Write composite averages to netCDF
avgtype = 1     ;0=Plot simple average; 1=Plot composite average WRT CPT height
cpzcomp = 17000 ;Used only if avgtype=1; Height at which to set composite CPT height
if(type.eq.1)then ;If we're writing to netCDF file, set cold point height to zero
   cpzcomp = 0
end if
minrad = 0      ;Minimum radius outside of which to pull sonde data
maxrad = 1000   ;Maximum radius within which to pull sonde data
dr     = 100    ;Width of radial bands
mincpz = 14000  ;Minimum cold-point tropopause height allowed (m)
ymin   = -5;14.
ymax   = 5.;21.
g = 9.80665
rads  = ispan(minrad,maxrad,dr)
nrads = (maxrad-minrad)/dr

if(intype.eq.1.and.type.eq.1)then
   print( "ERROR: intype=1 and type=1, which would overwrite input file.")
   print( "Set intype and type such that they're different.")
   exit
end if

if(machine.eq."albany")then
   if(intype.eq.0)
      infile = "/nfs/jm13/pduran/rawinsondes/sparc/rawinsondes_0-1000km.nc"
   else
      infile = "/nfs/jm13/pduran/rawinsondes/sparc/rawinsondes_"+minrad+"-"\
                +maxrad+"km_cpzcomposite_"+cat+".nc"
   end if
   outdir    = "/nfs/jm13/pduran/web/trop/"+cat+"/"
   ncoutfile = "/nfs/jm13/pduran/rawinsondes/sparc/rawinsondes_"+minrad+"-"\
                +maxrad+"km_cpzcomposite_"+cat+".nc"
else
   if(intype.eq.0)then
      infile = "/home/patrick/research/data/rawinsonde/rawinsondes_0-1000km.nc"
   else
      infile = "/home/patrick/research/data/rawinsonde/rawinsondes_"+minrad+"-"\
                +maxrad+"km_cpzcomposite_"+cat+".nc"
   end if
   outdir    = "/home/patrick/research/plots/trop/"+cat+"/"
   ncoutfile = "/home/patrick/research/data/rawinsonde/rawinsondes_"+minrad+"-"\
                +maxrad+"km_cpzcomposite_"+cat+".nc"
end if

system("mkdir "+outdir)

;Read input file
print( "Reading input file "+infile+"   "+systemfunc("date"))
ncfile = addfile(infile,"r")
if(intype.eq.0)then
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
   stormtypein= ncfile->stormtype
   stormlatin = ncfile->stormlat
   stormlonin = ncfile->stormlon
   vmaxin     = ncfile->vmax
   pminin     = ncfile->pmin
   irbtin     = ncfile->irbt
else
   n2_avg = ncfile->n2_avg
   tc_avg = ncfile->tc_avg
   p_avg  = ncfile->p_avg
   nsondes= ncfile->nsondes
   radius = ncfile->radius
end if

if(intype.eq.0)then
   ;Total number of sondes in the entire dataset
   totnumsondes = dimsizes(droptimein)
   
   ;Compute squared Brunt-Vaisala frequency
   tk = tcin+273.15
   th = pot_temp(presin,tk,-1,False)
   dthdz = center_finite_diff_n(th,heightin,False,0,1)  ;Units K/m
   ;A couple of sondes in the dataset have outlier tk values that
   ;cause dthdz to be nan. Replace these nans with missing values.
   replace_ieeenan(dthdz,dthdz@_FillValue,0)
   n2in = (g/th)*dthdz*10000 ;Units 10^-4 s^-2
   replace_ieeenan(n2in,n2in@_FillValue,0)

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
      tcshift = new((/100000,nhghts*2-1/),"float")
      tcshift!0 = "sondes in category"
      tcshift!1 = "comphght"
      tcshift&comphght = ispan(max(heightin)*-1,max(heightin),dz)*1.
      n2shift = tcshift
      pshift  = tcshift
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
   
   ;Initialize average arrays for compositing
   if(avgtype.eq.1)then
      tc_avg := new((/nrads,nhghts*2-1/),"float")
      n2_avg := new((/nrads,nhghts*2-1/),"float")
       p_avg := new((/nrads,nhghts*2-1/),"float")
   end if
   
   ;Get data only within the selected radial band
   do i=0,nrads-1
      unixi = i + 1
      minbinrad = rads(i)
      maxbinrad = rads(i+1)
      radinds  := ind(radiusin.ge.minbinrad.and.radiusin.lt.maxbinrad)
   
      droptime := droptimein(radinds) 
      droplat  := droplatin (radinds) 
      droplon  := droplonin (radinds) 
      height   := heightin 
      cpt      := cptin      (radinds) 
      cpz      := cpzin      (radinds) 
      tophght  := tophghtin  (radinds) 
      pres     := presin     (radinds,:) 
      tc       := tcin       (radinds,:) 
      sondelat := sondelatin (radinds,:) 
      sondelon := sondelonin (radinds,:) 
      radius   := radiusin   (radinds) 
      azimuth  := azimuthin  (radinds) 
      stormtype:= stormtypein(radinds) 
      stormlat := stormlatin (radinds) 
      stormlon := stormlonin (radinds) 
      vmax     := vmaxin     (radinds) 
      pmin     := pminin     (radinds) 
      n2       := n2in       (radinds,:)
      irbt     := irbtin     (radinds)

      n2!0 = "radius"
      n2!1 = "height"
      n2&height = height
   
      ;Get indices for each storm category
      if(intorir.eq.0)then
         if(cat.eq."td")then
            catinds := ind(vmax.le.33.and.stormtype.eq.0)
         else if(cat.eq."ts")then
            catinds := ind(vmax.ge.34.and.vmax.lt.64.and.stormtype.eq.0)
         else if(cat.eq."td+ts")then
            catinds := ind(vmax.lt.64.and.stormtype.eq.0)
         else if(cat.eq."hur")then
            catinds := ind(vmax.ge.64.and.stormtype.eq.0)
         else if(cat.eq."cat1+2")then
            catinds := ind(vmax.ge.64.and.vmax.le.95.and.stormtype.eq.0)
         else if(cat.eq."maj")then
            catinds := ind(vmax.gt.95.and.stormtype.eq.0)
         else if(cat.eq."all")then
            catinds := ind(vmax.ge.0.and.stormtype.eq.0)
         end if
         end if
         end if
         end if
         end if
         end if
         end if
      else
         if(irscene.eq."coldest")then
            catinds := ind(irbt.le.coldcutoff)
         else if(irscene.eq."cold")then
            catinds := ind(irbt.gt.coldcutoff.and.irbt.le.warmcutoff)
         else if(irscene.eq."warm")then
            catinds := ind(irbt.gt.warmcutoff)
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
      else if(avgtype.eq.1) ;Compute composite average with respect to cold-point height
         nsondesincat = dimsizes(catinds)
         do j=0,nsondesincat-1    ;Loop through all sondes in this intensity category
            unixj = j + 1
            print( "Processing sonde "+unixj+" of "+nsondesincat+" for radbin "+unixi+\
                   " of "+nrads+"   "+systemfunc("date"))
            catind = catinds(j)
            cpzcur = cpz_cat(i,j)      ;Current cold point height
            if(.not.ismissing(cpzcur))
               maxhghtshift = max(height)-cpzcur
               minhghtshift = -cpzcur+dz
               tccur = tc  (catind,:)    ;T profile for this sonde
               n2cur = n2  (catind,:)    ;N^2 profile for this sonde
                pcur = pres(catind,:)    ;P profile for this sonde
               tcshift(j,{0:maxhghtshift}) = (/tccur({cpzcur:})/)
               tcshift(j,{minhghtshift:0}) = (/tccur({:cpzcur-dz})/)
               n2shift(j,{0:maxhghtshift}) = (/n2cur({cpzcur:})/)
               n2shift(j,{minhghtshift:0}) = (/n2cur({:cpzcur-dz})/)
                pshift(j,{0:maxhghtshift}) = (/ pcur({cpzcur:})/)
                pshift(j,{minhghtshift:0}) = (/ pcur({:cpzcur-dz})/)
            end if
         end do
         tc_avg(i,:) = dim_avg_n_Wrap(tcshift,0)
         n2_avg(i,:) = dim_avg_n_Wrap(n2shift,0)
          p_avg(i,:) = dim_avg_n_Wrap( pshift,0)
      end if
      end if
      ;Determine cold point tropopause temperature and height using mean fields
      if(avgtype.eq.0)then
         cpz_avg(i) = height(minind(tc_avg(i,:)))/1000. ;Convert to km
         cpt_avg(i) = min(tc_avg(i,:))
      else
         cpz_avg(i) = cpzcomp/1000.
      end if
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
   if(avgtype.eq.0)then
      n2_avg!1 = "height"
      n2_avg&radius = midpoints
      n2_avg&height = heightkm
   else
      n2_avg&radius = midpoints
      n2_avg&comphght = (n2_avg&comphght+cpzcomp)/1000.
   end if
   copy_VarCoords(n2_avg,tc_avg)
   copy_VarCoords(n2_avg, p_avg)
end if

if(intype.eq.0.and.type.eq.0)then
   ;Compute and plot probability distributions
   cpzopt = True
   cpzopt@bin_min = 14.
   cpzopt@bin_max = 19.
   cpzopt@bin_nice= True
   pdf_cpz = pdfx(cpz_cat,21,cpzopt)
   cptopt = cpzopt
   cptopt@bin_min = -90
   cptopt@bin_max = -50
   pdf_cpt = pdfx(cpt_cat,61,cptopt)

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
end if

;If desired, make some plots
if(type.eq.0)then
   ;-------------------;
   ;Cross-section plots;
   ;-------------------;
   out_suffix = minrad+"-"+maxrad+"km_"+ymin+"-"+ymax+"km"
   ;Transpose arrays so that height is leftmost dimension
   n2_avg_plot = transpose(n2_avg)
   tc_avg_plot = transpose(tc_avg)

   ;Set output file names
   if(avgtype.eq.0)then
      outfilen2 = outdir+"stab+cp_trop_hght_"+out_suffix
      outfiletc = outdir+"temp+cp_trop_hght_"+out_suffix
   else
      outfilen2 = outdir+"stab+cp_trop_hght_composite"+out_suffix
      outfiletc = outdir+"temp+cp_trop_hght_composite"+out_suffix
   end if

   ;STATIC STABILITY AND COLD-POINT TROPOPAUSE HEIGHT
   wks_stab = gsn_open_wks("png",outfilen2)
   gsn_define_colormap(wks_stab,"WhiteGreen")
   stabres = True
   stabres@cnFillOn = True
   stabres@cnLevelSelectionMode = "ManualLevels"
   stabres@cnMinLevelValF = 0
   stabres@cnMaxLevelValF = 10.0
   stabres@cnLevelSpacingF = 1.0
   stabres@cnLineLabelsOn = False
   stabres@cnLinesOn = False
   stabres@gsnLeftStringFontHeightF = 0.014
   stabres@gsnRightStringFontHeightF = stabres@gsnLeftStringFontHeightF
   stabres@gsnLeftString = "N~S~2~N~ (10~S~-4~N~ s~S~-2~N~) and cold point height"
   stabres@gsnRightString = "Category: "+str_upper(cat)+"        "+minrad+"-"+maxrad+" km"
   stabres@tiXAxisString = "Radius (km)"
   if(avgtype.eq.0)then
      stabres@tiYAxisString = "Height (km)"
   else
      stabres@tiYAxisString = "Height relative to cold point (km)"
   end if
   stabres@lbLabelBarOn = True
   stabres@lbAutoManage = False
   stabres@lbOrientation = "Vertical"
   stabres@cnInfoLabelOn = False
   ;stabres@trYMinF = 9
   ;stabres@trYMaxF = 25
   stabres@trYMinF = ymin
   stabres@trYMaxF = ymax
   stabres@tmYROn = False
   stabres@tmXTOn = False
   stabres@gsnMaximize = True
   stabres@gsnDraw = False
   stabres@gsnFrame = False
   stab_img = gsn_csm_contour(wks_stab,n2_avg_plot,stabres)
   ;Add number of sondes along bottom axis
   txres=True
   txres@txFontHeightF = 0.02
   txx = radius
   txy = ymin+0.5
   text= gsn_add_text(wks_stab,stab_img,tostring(nsondes),txx,txy,txres)
   tropres = True
   tropres@gsLineColor = "yellow"
   tropres@gsLineThicknessF = 18
   ;trophghts = gsn_add_polyline(wks_stab,stab_img,cpz_avg&radius,cpz_avg,tropres)
   print( "Writing N^2 image to "+outfilen2+".png")
   draw(stab_img)
   frame(wks_stab)

   ;TEMPERATURE
   wks_tc = gsn_open_wks("png",outfiletc)
   gsn_define_colormap(wks_tc,"rainbow")
   tcres = stabres
   tcres@gsnLeftString = "Temperature ~S~o~N~C"
   tcres@cnMinLevelValF = -80
   tcres@cnMaxLevelValF = -40
   tcres@cnLevelSpacingF = 1.
   tc_img = gsn_csm_contour(wks_tc,tc_avg_plot,tcres)
   ;Add number of sondes along bottom axis
   text= gsn_add_text(wks_tc,tc_img,tostring(nsondes),txx,txy,txres)
   ;Overlay contours for coldest brightness temperatures
   tccontres = stabres
   tccontres@gsnLeftString = " "
   tccontres@cnFillOn = False
   tccontres@cnMinLevelValF = -90
   tccontres@cnMaxLevelValF = -74
   tccontres@cnLevelSpacingF = 1 
;   tccontres@cnLineLabelsOn = True
   tccontres@cnLinesOn = True
   tccontres@cnLineThicknessF = 4.
   tccontres@cnLineColor = "yellow"
   tccontimg = gsn_csm_contour(wks_tc,tc_avg_plot,tccontres)
   overlay(tc_img,tccontimg)
   print( "Writing temperature image to "+outfilen2+".png")
   draw(tc_img)
   frame(wks_tc)
end if

;Write to composite output file, if desired
if(intype.eq.0.and.type.eq.1)
   ;Set attributes for each variable
   comprad  = n2_avg&radius
   comphght = n2_avg&comphght
   comprad@long_name = "distance from TC center"
   comprad@units = "km"
   comphght@long_name = "height above or below cold-point tropopause"
   comphght@units = "km"
   comphght@_FillValue = n2_avg@_FillValue
   nsondes@long_name = "Number of rawinsondes included in the average (must have data up to the cold point)"
   n2_avg@long_name = "squared Brunt-Vaisala frequency"
   n2_avg@units = "10^-4 s^-2"
   tc_avg@long_name = "temperature"
   tc_avg@units = "deg C"
    p_avg@long_name = "pressure"
    p_avg@units = "hPa"

   print( "Writing variables to output file...   "+systemfunc("date"))
   print(" "+ncoutfile)
   system("rm -f "+ncoutfile)
   fout = addfile(ncoutfile,"c")
   setfileoption(fout,"DefineMode",True)

   ;Define global attributes
   if(intorir.eq.0)then
      if(cat.eq.td)then
         cgroup = "tropical depressions"
      else if(cat.eq."ts")then
         cgroup = "tropical storms"
      else if(cat.eq."td+ts")then
         cgroup = "tropical depressions and tropical storms"
      else if(cat.eq."hur")then
         cgroup = "all hurricanes"
      else if(cat.eq."cat1+2")then
         cgroup = "category 1 and 2 hurricanes"
      else if(cat.eq."maj")then
         cgroup = "major hurricanes"
      else if(cat.eq."all")then
         cgroup = "all sondes"
      end if
      end if
      end if
      end if
      end if
      end if
      end if
   else
      if(irscene.eq."coldest")then
         cgroup = "brightness temperatures colder than "+coldcutoff
      else if(irscene.eq."cold")then
         cgroup = "brightness temperatures between "+warmcutoff+" and "+coldcutuff
      else if(irscene.eq."warm")then
         cgroup = "brightness temperatures warmer than "+warmcutoff
      end if
      end if
      end if
   end if
   fileatts = True
   fileatts@minrad = minrad+" km"
   fileatts@maxrad = maxrad+" km"
   fileatts@dr     = dr+" km"
   fileatts@mincpz = mincpz+" m"
   fileatts@description = "Contains composite averages computed with respect to"+\
			  " the cold-point tropopause height. comphght=0 is the"+\
                          " value at the cold point. These averages are"+\
                          " computed for "+cgroup
   fileattdef(fout,fileatts)

   ;Define variable attributes
   dimNames = (/"radius","comphght"/)
   dimSizes = (/dimsizes(comprad),dimsizes(comphght)/)
   dimUnlim = (/False,False/)
   filedimdef(fout,dimNames,dimSizes,dimUnlim)

   ;Define coordinate variables for writing
   filevardef(fout,"comphght",typeof(comphght),getvardims(comphght))
   filevardef(fout,"radius",typeof(comprad),getvardims(comprad))
   filevardef(fout,"nsondes",typeof(nsondes),getvardims(nsondes))
   filevarattdef(fout,"comphght",comphght)
   filevarattdef(fout,"radius",comprad)
   filevarattdef(fout,"nsondes",comprad)

   ;Define the variables to write
   var_names = (/"n2_avg","tc_avg","p_avg"/)
   var_types = (/"float","float","float"/)
   filevardef(fout,var_names,var_types,dimNames)
   filevarattdef(fout,"n2_avg",n2_avg)
   filevarattdef(fout,"tc_avg",tc_avg)
   filevarattdef(fout,"p_avg" , p_avg)

   ;Write variables to output file
   ;NOTE: The (/xxx/) syntax here tells NCL to output only the data
   ;values to the predefined locations in the file. This greatly
   ;inproves the efficiency of the writing routine.
   printVarSummary( comphght)
   fout->comphght = (/comphght/)
   fout->radius = (/comprad/)
   fout->nsondes= (/nsondes/)
   fout->n2_avg = (/n2_avg/)
   fout->tc_avg = (/tc_avg/)
   fout->p_avg  = (/ p_avg/)
end if
