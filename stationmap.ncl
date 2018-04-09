;---------------------------------------------
; stationmap.ncl
; Patrick Duran
; 09 April 2018
; 
; Reads in ListOfStations.dat, gets the number
; of rawinsondes in the hurricane dataset from
; rawinsondes_0-1000km for each station, and
; plots the numbers on a map at the location
; of each station.
;---------------------------------------------

machine = "personal"  ;Set to either "personal" or "albany"
if(machine.eq."personal")then
   stnlist   = "/home/patrick/research/data/rawinsonde/ListOfStations.dat"
   sondelist = "/home/patrick/research/data/rawinsonde/rawinsondes_0-1000km"
   outdir    = "/home/patrick/research/plots/"
else
   stnlist   = "/jm13/pduran/hurricanesondes/ListOfStations.dat"
   sondelist = "/jm13/pduran/hurricanesondes/rawinsondes_0-1000km"
   outdir    = "/jm13/pduran/web/trop/"
end if

stids = systemfunc("cat "+stnlist+" | sed '1d' | awk -F',' '{print $1}'")
lats  = tofloat(systemfunc("cat "+stnlist+" | sed '1d' | awk -F',' '{print $6}'"))
lons  = tofloat(systemfunc("cat "+stnlist+" | sed '1d' | awk -F',' '{print $7}'"))

;Loop through all of the stations and get number
;of sondes from each station included in the dataset.
nstns = dimsizes(stids)
nsondes = new(nstns,"string")
do i=0,nstns-1
   unixi = i + 1
   print( "Processing station "+unixi+" of "+nstns+"   "+systemfunc("date"))
   stn = stids(i)
   nsondes(i) = systemfunc("grep -i "+stn+"- "+sondelist+" | wc | awk '{print $1}'")
end do

;Plot the map
wkstype = "png"
wks = gsn_open_wks(wkstype,outdir+"stationmap")
res = True
res@gsnMaximize = True
res@mpFillOn = False
res@mpLimitMode = "LatLon"
res@mpMinLatF = 10
res@mpMaxLatF = 50
res@mpMinLonF = -110
res@mpMaxLonF = -55
res@mpOutlineBoundarySets = "GeophysicalAndUSStates"
res@tmYROn = False
res@tmXTOn = False
res@gsnLeftString = "Number of sondes by station location"
res@gsnLeftStringOrthogonalPosF = -0.003
res@tiXAxisString = "Degrees Longitude"
res@tiYAxisString = "Degrees Latitude"
res@tiXAxisFontHeightF = 0.02
res@tiYAxisFontHeightF = res@tiXAxisFontHeightF
res@tiXAxisOffsetYF = -0.05
res@tiYAxisOffsetXF = -0.07
res@gsnDraw = False
res@gsnFrame= False
map = gsn_csm_map(wks,res)

;Overlay text to indicate number of sondes at each station
txres = True
txres@txFont = 1
txres@txFontColor = "blue"
txres@txFontHeightF = 0.01
txres@txFontThicknessF = 4.
txres@txJust = "centercenter"
tx = gsn_add_text(wks,map,nsondes,lons,lats,txres)
draw(map)
frame(wks)
