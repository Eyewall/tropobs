;---------------------------------------------
; stormtracks.ncl
; Patrick Duran
; 09 April 2018
; 
; Reads in rawinsondes_0-1000km.nc, gets storm
; lats, lons, and vmax, and plots color-coded
; dots to indicate the storm center positions
; and maximum intensity for each sonde obs.
;---------------------------------------------

machine = "albany"  ;Set to either "personal" or "albany"
if(machine.eq."personal")then
   infile   = "/home/patrick/research/data/rawinsonde/rawinsondes_0-1000km.nc"
   outdir    = "/home/patrick/research/plots/"
else
   infile   = "/jm13/pduran/rawinsondes/sparc/rawinsondes_0-1000km.nc"
   outdir    = "/jm13/pduran/web/trop/"
end if

;Read variables from input file
ncfile = addfile(infile,"r")
lats = ncfile->stormlat
lons = ncfile->stormlon
vmax = ncfile->vmax
stormtype = ncfile->stormtype

;Get indices for each storm category
tdinds = ind(vmax.lt.34.and.stormtype.eq.0)
tsinds = ind(vmax.ge.34.and.vmax.lt.64.and.stormtype.eq.0)
hurinds= ind(vmax.ge.64.and.stormtype.eq.0)

;Plot the map
wkstype = "png"
wks = gsn_open_wks(wkstype,outdir+"stormtracks")
res = True
res@gsnMaximize = True
res@mpFillOn = False
res@mpLimitMode = "LatLon"
res@mpMinLatF = 10
res@mpMaxLatF = 50
res@mpMinLonF = -110
res@mpMaxLonF = -50
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
res@tmYROn = False
res@tmXTOn = False
res@gsnLeftString = "Storm center positions when observed by rawinsondes"
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = -0.003
res@tiXAxisString = "Degrees Longitude"
res@tiYAxisString = "Degrees Latitude"
res@tiXAxisFontHeightF = 0.015
res@tiYAxisFontHeightF = res@tiXAxisFontHeightF
res@tiXAxisOffsetYF = -0.05
res@tiYAxisOffsetXF = -0.07
res@gsnDraw = False
res@gsnFrame= False
map = gsn_csm_map(wks,res)

;Overlay dots representing storm center locations
dotres = True
dotres@gsMarkerColor = "blue"
tddots  = gsn_add_polymarker(wks,map,lons(tdinds),lats(tdinds),dotres)
dotres@gsMarkerColor = "orange"
tsdots  = gsn_add_polymarker(wks,map,lons(tsinds),lats(tsinds),dotres)
dotres@gsMarkerColor = "red"
hurdots = gsn_add_polymarker(wks,map,lons(hurinds),lats(hurinds),dotres)
draw(map)
frame(wks)
