  remark "geqdsk to ifile"

  read "d3.bas"
  read "sewall.bas"

  character*32 gname = "eqdsk" 
  character*32 fileid = "shot_time"
  character*32 inv_save_file = "shot_time_inv.sav"

  d3("eqdsk",0,0.995,1,501,0.001)

  fileid = trim(shotName) // "_" // format(shotTime*1e3,-4,0,1) 
  inv_save_file = trim(fileid) // "_inv.sav"

  restore ^inv_save_file
  call teq_inv
  integer mymap=257		#number of poloidal grid points
  integer mynht=600		#max number of equilibrium solve interations
  eq.map=mymap
  eq.nht=mynht
  eq.map;eq.nht
  generate

  weqdsk("i","_test")

  quit(1)
