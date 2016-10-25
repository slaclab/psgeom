; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; Manually optimized with hdfsee
; This geometry started with a geometry that workde for indexing 2009 PIS data 
; collected at AMO. 

adu_per_eV = 0.078
photon_energy = /LCLS/photon_energy_eV
clen = /LCLS/detector_1/EncoderValue
;clen = 143.0 for runs 15 - 25

data = /entry_1/instrument_1/detector_1/detector_corrected/data
mask = /entry_1/instrument_1/detector_1/detector_corrected/mask
mask_good = 0x0000
mask_bad = 0xffff
dim0 = %
dim1 = ss
dim2 = fs

rigid_group_q0 = top
rigid_group_q1 = bottom

rigid_group_a0 = top
rigid_group_a1 = bottom

rigid_group_collection_quadrants = q0,q1
rigid_group_collection_asics = a0,a1

; These default values will be used unless overridden by the per-panel values

top/min_fs = 512
top/max_fs = 1023
top/min_ss = 0
top/max_ss = 1023
top/corner_x = -519
top/corner_y = 6
top/fs = +0.000000x +1.000000y
top/ss = +1.000000x +0.000000y
top/res = 13333.3  ; 75 micron pixel size
top/badrow_direction = -
top/coffset = 0.47006

bottom/min_fs = 0
bottom/max_fs = 511
bottom/min_ss = 0
bottom/max_ss = 1023
bottom/corner_x = -515
bottom/corner_y = -646
bottom/fs = +0.000000x +1.000000y
bottom/ss = +1.000000x +0.000000y
bottom/res = 13333.3  ; 75 micron pixel size
bottom/badrow_direction = -
bottom/coffset = 0.4725

