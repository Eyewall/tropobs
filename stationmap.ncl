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

machine = "albany"  ;Set to either "personal" or "albany"
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

;Split input strings to get only the storm types and station numbers
nlines  = stringtoint(systemfunc("wc "+sondelist+" | awk '{print $1}'"))
inlines = asciiread(sondelist,nlines,"string")
linesplit = str_split_csv(inlines," ",3)
stormtypes = stringtoint(linesplit(:,6))
dirsplit = str_split_csv(linesplit(:,0),"/",3)
filenamesplit = str_split_csv(dirsplit(:,7),"-",3)
stnnums = stringtoint(filenamesplit(:,0))

;Loop through all of the stations and get number
;of sondes from each station included in the dataset.
nstns = dimsizes(stids)
nsondes = new(nstns,"string")
do i=0,nstns-1
   unixi = i + 1
   print( "Processing station "+unixi+" of "+nstns+"   "+systemfunc("date"))
   stn = stringtoint(stids(i))
   ;None of the observations from these three stations are in tropical systems,
   ;but for some reason that I haven't been able to figure out, the script
   ;still wants to count them. Set all to a space so that the numbers are not
   ;plotted. Else, set each station to the appropriate number.
   if(stn.eq.14918.or.stn.eq.23062.or.stn.eq.23066)then
      nsondes(i) = " "
   else
      nsondes(i) = dimsizes(ind(stnnums.eq.stn.and.stormtypes.eq.0))
   end if
end do
totnumsondes = sum(toint(nsondes))

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
txres@txFont = 22
txres@txFontColor = "blue"
txres@txFontHeightF = 0.015
txres@txFontThicknessF = 4.
txres@txJust = "centercenter"
tx = gsn_add_text(wks,map,nsondes,lons,lats,txres)

;Add total number of sondes
totnumtext = "Total: "+totnumsondes
txres@txFont = 22
txres@txFontHeightF = 0.03
;txres@txFontThicknessF = 10.
txres@txJust = "BottomLeft"
tottx = gsn_add_text(wks,map,totnumtext,res@mpMinLonF+.3,res@mpMinLatF+.3,txres)
draw(map)
frame(wks)
