; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Optimized panel offsets can be found at the end of the file
adu_per_eV = 0.0042  ; 42 adu/keV for highest gain mode
clen =  0.1015
photon_energy = 15000
res = 13333.3	; 75 um pixels

data = /data/data


bad_s0/min_fs = 0
bad_s0/max_fs = 0
bad_s0/min_ss = 0
bad_s0/max_ss = 1023

bad_s1/min_fs = 255
bad_s1/max_fs = 256
bad_s1/min_ss = 0
bad_s1/max_ss = 1023

bad_s2/min_fs = 511
bad_s2/max_fs = 512
bad_s2/min_ss = 0
bad_s2/max_ss = 1023

bad_s3/min_fs = 767
bad_s3/max_fs = 768
bad_s3/min_ss = 0
bad_s3/max_ss = 1023

bad_s4/min_fs = 1023
bad_s4/max_fs = 1023
bad_s4/min_ss = 0
bad_s4/max_ss = 1023

bad_m0/min_fs = 0
bad_m0/max_fs = 1023
bad_m0/min_ss = 0
bad_m0/max_ss = 0

bad_m1/min_fs = 0
bad_m1/max_fs = 1023
bad_m1/min_ss = 255
bad_m1/max_ss = 256

bad_m2/min_fs = 0
bad_m2/max_fs = 1023
bad_m2/min_ss = 511
bad_m2/max_ss = 512

bad_m3/min_fs = 0
bad_m3/max_fs = 1023
bad_m3/min_ss = 767
bad_m3/max_ss = 768

bad_m4/min_fs = 0
bad_m4/max_fs = 1023
bad_m4/min_ss = 1023
bad_m4/max_ss = 1023

rigid_group_p1 = p1a1,p1a2,p1a3,p1a4,p1a5,p1a6,p1a7,p1a8
rigid_group_p2 = p2a1,p2a2,p2a3,p2a4,p2a5,p2a6,p2a7,p2a8
rigid_group_collection_det = p1,p2

p1a1/corner_x = -964
p1a1/corner_y = 31
p1a1/fs = +1.000000x +0.000000y
p1a1/ss = +0.000000x +1.000000y
p1a1/min_fs = 0
p1a1/max_fs = 255
p1a1/min_ss = 0
p1a1/max_ss = 255

p1a2/corner_x = -706
p1a2/corner_y = 31
p1a2/fs = +1.000000x +0.000000y
p1a2/ss = +0.000000x +1.000000y
p1a2/min_fs = 256
p1a2/max_fs = 511
p1a2/min_ss = 0
p1a2/max_ss = 255

p1a3/corner_x = -448
p1a3/corner_y = 31
p1a3/fs = +1.000000x +0.000000y
p1a3/ss = +0.000000x +1.000000y
p1a3/min_fs = 512
p1a3/max_fs = 767
p1a3/min_ss = 0
p1a3/max_ss = 255

p1a4/corner_x = -190
p1a4/corner_y = 31
p1a4/fs = +1.000000x +0.000000y
p1a4/ss = +0.000000x +1.000000y
p1a4/min_fs = 768
p1a4/max_fs = 1023
p1a4/min_ss = 0
p1a4/max_ss = 255

p1a5/corner_x = -964
p1a5/corner_y = 289
p1a5/fs = +1.000000x +0.000000y
p1a5/ss = +0.000000x +1.000000y
p1a5/min_fs = 0
p1a5/max_fs = 255
p1a5/min_ss = 256
p1a5/max_ss = 511

p1a6/corner_x = -706
p1a6/corner_y = 289
p1a6/fs = +1.000000x +0.000000y
p1a6/ss = +0.000000x +1.000000y
p1a6/min_fs = 256
p1a6/max_fs = 511
p1a6/min_ss = 256
p1a6/max_ss = 511

p1a7/corner_x = -448
p1a7/corner_y = 289
p1a7/fs = +1.000000x +0.000000y
p1a7/ss = +0.000000x +1.000000y
p1a7/min_fs = 512
p1a7/max_fs = 767
p1a7/min_ss = 256
p1a7/max_ss = 511

p1a8/corner_x = -190
p1a8/corner_y = 289
p1a8/fs = +1.000000x +0.000000y
p1a8/ss = +0.000000x +1.000000y
p1a8/min_fs = 768
p1a8/max_fs = 1023
p1a8/min_ss = 256
p1a8/max_ss = 511

p2a1/corner_x = -964
p2a1/corner_y = -519
p2a1/fs = +1.000000x +0.000000y
p2a1/ss = +0.000000x +1.000000y
p2a1/min_fs = 0
p2a1/max_fs = 255
p2a1/min_ss = 512
p2a1/max_ss = 767

p2a2/corner_x = -706
p2a2/corner_y = -519
p2a2/fs = +1.000000x +0.000000y
p2a2/ss = +0.000000x +1.000000y
p2a2/min_fs = 256
p2a2/max_fs = 511
p2a2/min_ss = 512
p2a2/max_ss = 767

p2a3/corner_x = -448
p2a3/corner_y = -519
p2a3/fs = +1.000000x +0.000000y
p2a3/ss = +0.000000x +1.000000y
p2a3/min_fs = 512
p2a3/max_fs = 767
p2a3/min_ss = 512
p2a3/max_ss = 767

p2a4/corner_x = -190
p2a4/corner_y = -519
p2a4/fs = +1.000000x +0.000000y
p2a4/ss = +0.000000x +1.000000y
p2a4/min_fs = 768
p2a4/max_fs = 1023
p2a4/min_ss = 512
p2a4/max_ss = 767

p2a5/corner_x = -964
p2a5/corner_y = -261
p2a5/fs = +1.000000x +0.000000y
p2a5/ss = +0.000000x +1.000000y
p2a5/min_fs = 0
p2a5/max_fs = 255
p2a5/min_ss = 768
p2a5/max_ss = 1023

p2a6/corner_x = -706
p2a6/corner_y = -261
p2a6/fs = +1.000000x +0.000000y
p2a6/ss = +0.000000x +1.000000y
p2a6/min_fs = 256
p2a6/max_fs = 511
p2a6/min_ss = 768
p2a6/max_ss = 1023

p2a7/corner_x = -448
p2a7/corner_y = -261
p2a7/fs = +1.000000x +0.000000y
p2a7/ss = +0.000000x +1.000000y
p2a7/min_fs = 512
p2a7/max_fs = 767
p2a7/min_ss = 768
p2a7/max_ss = 1023

p2a8/corner_x = -190
p2a8/corner_y = -261
p2a8/fs = +1.000000x +0.000000y
p2a8/ss = +0.000000x +1.000000y
p2a8/min_fs = 768
p2a8/max_fs = 1023
p2a8/min_ss = 768
p2a8/max_ss = 1023
