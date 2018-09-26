; geoptimiser -i r0020_0.stream -g test4_opt5.geom -o test4_opt6.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4_opt4.geom -o test4_opt5.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4_opt3.geom -o test4_opt4.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4_opt2.geom -o test4_opt3.geom -c quadrants -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4_opt1.geom -o test4_opt2.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4_opt.geom -o test4_opt1.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; geoptimiser -i r0020_0.stream -g test4.geom -o test4_opt.geom -c asics -q quadrants
; Optimized panel offsets can be found at the end of the file
; Manually optimized with hdfsee
;Automatically generated from calibration data

clen =  /LCLS/detector0-EncoderValue
coffset = 0.50
photon_energy = /LCLS/photon_energy_eV
res = 9097.52
adu_per_eV = 0.00338
;data = /data/peakpowder

; The following lines define "rigid groups" which express the physical
; construction of the detector.  This is used when refining the detector
; geometry.

data = /entry_1/data_1/data
dim0 = %
dim1 = ss
dim2 = fs

rigid_group_q0 = q0a0,q0a1,q0a2,q0a3,q0a4,q0a5,q0a6,q0a7,q0a8,q0a9,q0a10,q0a11,q0a12,q0a13,q0a14,q0a15
rigid_group_q1 = q1a0,q1a1,q1a2,q1a3,q1a4,q1a5,q1a6,q1a7,q1a8,q1a9,q1a10,q1a11,q1a12,q1a13,q1a14,q1a15
rigid_group_q2 = q2a0,q2a1,q2a2,q2a3,q2a4,q2a5,q2a6,q2a7,q2a8,q2a9,q2a10,q2a11,q2a12,q2a13,q2a14,q2a15
rigid_group_q3 = q3a0,q3a1,q3a2,q3a3,q3a4,q3a5,q3a6,q3a7,q3a8,q3a9,q3a10,q3a11,q3a12,q3a13,q3a14,q3a15

rigid_group_a0 = q0a0,q0a1
rigid_group_a1 = q0a2,q0a3
rigid_group_a2 = q0a4,q0a5
rigid_group_a3 = q0a6,q0a7
rigid_group_a4 = q0a8,q0a9
rigid_group_a5 = q0a10,q0a11
rigid_group_a6 = q0a12,q0a13
rigid_group_a7 = q0a14,q0a15
rigid_group_a8 = q1a0,q1a1
rigid_group_a9 = q1a2,q1a3
rigid_group_a10 = q1a4,q1a5
rigid_group_a11 = q1a6,q1a7
rigid_group_a12 = q1a8,q1a9
rigid_group_a13 = q1a10,q1a11
rigid_group_a14 = q1a12,q1a13
rigid_group_a15 = q1a14,q1a15
rigid_group_a16 = q2a0,q2a1
rigid_group_a17 = q2a2,q2a3
rigid_group_a18 = q2a4,q2a5
rigid_group_a19 = q2a6,q2a7
rigid_group_a20 = q2a8,q2a9
rigid_group_a21 = q2a10,q2a11
rigid_group_a22 = q2a12,q2a13
rigid_group_a23 = q2a14,q2a15
rigid_group_a24 = q3a0,q3a1
rigid_group_a25 = q3a2,q3a3
rigid_group_a26 = q3a4,q3a5
rigid_group_a27 = q3a6,q3a7
rigid_group_a28 = q3a8,q3a9
rigid_group_a29 = q3a10,q3a11
rigid_group_a30 = q3a12,q3a13
rigid_group_a31 = q3a14,q3a15

rigid_group_collection_quadrants = q0,q1,q2,q3
rigid_group_collection_asics = a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31



q0a0/min_fs = 0
q0a0/min_ss = 0
q0a0/max_fs = 193
q0a0/max_ss = 184
q0a0/fs = -0.002384x +0.999997y
q0a0/ss = -0.999997x -0.002384y
q0a0/corner_x = 449.861
q0a0/corner_y = -32.5139
q0a0/no_index = 0

q0a1/min_fs = 194
q0a1/min_ss = 0
q0a1/max_fs = 387
q0a1/max_ss = 184
q0a1/fs = -0.002384x +0.999997y
q0a1/ss = -0.999997x -0.002384y
q0a1/corner_x = 449.391
q0a1/corner_y = 164.486
q0a1/no_index = 0

q0a2/min_fs = 0
q0a2/min_ss = 185
q0a2/max_fs = 193
q0a2/max_ss = 369
q0a2/fs = -0.004088x +0.999992y
q0a2/ss = -0.999992x -0.004088y
q0a2/corner_x = 243.787
q0a2/corner_y = -33.1614
q0a2/no_index = 0

q0a3/min_fs = 194
q0a3/min_ss = 185
q0a3/max_fs = 387
q0a3/max_ss = 369
q0a3/fs = -0.004088x +0.999992y
q0a3/ss = -0.999992x -0.004088y
q0a3/corner_x = 242.983
q0a3/corner_y = 163.839
q0a3/no_index = 0

q0a4/min_fs = 0
q0a4/min_ss = 370
q0a4/max_fs = 193
q0a4/max_ss = 554
q0a4/fs = -0.999993x -0.003915y
q0a4/ss = +0.003915x -0.999993y
q0a4/corner_x = 871.379
q0a4/corner_y = 363.347
q0a4/no_index = 0

q0a5/min_fs = 194
q0a5/min_ss = 370
q0a5/max_fs = 387
q0a5/max_ss = 554
q0a5/fs = -0.999993x -0.003915y
q0a5/ss = +0.003915x -0.999993y
q0a5/corner_x = 674.378
q0a5/corner_y = 362.576
q0a5/no_index = 0

q0a6/min_fs = 0
q0a6/min_ss = 555
q0a6/max_fs = 193
q0a6/max_ss = 739
q0a6/fs = -0.999990x -0.004527y
q0a6/ss = +0.004527x -0.999990y
q0a6/corner_x = 872.175
q0a6/corner_y = 157.181
q0a6/no_index = 0

q0a7/min_fs = 194
q0a7/min_ss = 555
q0a7/max_fs = 387
q0a7/max_ss = 739
q0a7/fs = -0.999990x -0.004527y
q0a7/ss = +0.004527x -0.999990y
q0a7/corner_x = 675.176
q0a7/corner_y = 156.288
q0a7/no_index = 0

q0a8/min_fs = 0
q0a8/min_ss = 740
q0a8/max_fs = 193
q0a8/max_ss = 924
q0a8/fs = +0.004078x -0.999992y
q0a8/ss = +0.999992x +0.004078y
q0a8/corner_x = 478.419
q0a8/corner_y = 790.079
q0a8/no_index = 0

q0a9/min_fs = 194
q0a9/min_ss = 740
q0a9/max_fs = 387
q0a9/max_ss = 924
q0a9/fs = +0.004078x -0.999992y
q0a9/ss = +0.999992x +0.004078y
q0a9/corner_x = 479.224
q0a9/corner_y = 593.079
q0a9/no_index = 0

q0a10/min_fs = 0
q0a10/min_ss = 925
q0a10/max_fs = 193
q0a10/max_ss = 1109
q0a10/fs = +0.007235x -0.999973y
q0a10/ss = +0.999973x +0.007235y
q0a10/corner_x = 688.567
q0a10/corner_y = 790.567
q0a10/no_index = 0

q0a11/min_fs = 194
q0a11/min_ss = 925
q0a11/max_fs = 387
q0a11/max_ss = 1109
q0a11/fs = +0.007235x -0.999973y
q0a11/ss = +0.999973x +0.007235y
q0a11/corner_x = 689.992
q0a11/corner_y = 593.571
q0a11/no_index = 0

q0a12/min_fs = 0
q0a12/min_ss = 1110
q0a12/max_fs = 193
q0a12/max_ss = 1294
q0a12/fs = -0.999970x -0.007773y
q0a12/ss = +0.007773x -0.999970y
q0a12/corner_x = 446.761
q0a12/corner_y = 765.079
q0a12/no_index = 0

q0a13/min_fs = 194
q0a13/min_ss = 1110
q0a13/max_fs = 387
q0a13/max_ss = 1294
q0a13/fs = -0.999970x -0.007773y
q0a13/ss = +0.007773x -0.999970y
q0a13/corner_x = 249.766
q0a13/corner_y = 763.55
q0a13/no_index = 0

q0a14/min_fs = 0
q0a14/min_ss = 1295
q0a14/max_fs = 193
q0a14/max_ss = 1479
q0a14/fs = -0.999991x -0.004068y
q0a14/ss = +0.004068x -0.999991y
q0a14/corner_x = 445.988
q0a14/corner_y = 556.557
q0a14/no_index = 0

q0a15/min_fs = 194
q0a15/min_ss = 1295
q0a15/max_fs = 387
q0a15/max_ss = 1479
q0a15/fs = -0.999991x -0.004068y
q0a15/ss = +0.004068x -0.999991y
q0a15/corner_x = 248.991
q0a15/corner_y = 555.756
q0a15/no_index = 0

q1a0/min_fs = 388
q1a0/min_ss = 0
q1a0/max_fs = 581
q1a0/max_ss = 184
q1a0/fs = -0.999997x -0.002445y
q1a0/ss = +0.002445x -0.999997y
q1a0/corner_x = 41.6176
q1a0/corner_y = 438.051
q1a0/no_index = 0

q1a1/min_fs = 582
q1a1/min_ss = 0
q1a1/max_fs = 775
q1a1/max_ss = 184
q1a1/fs = -0.999997x -0.002445y
q1a1/ss = +0.002445x -0.999997y
q1a1/corner_x = -155.381
q1a1/corner_y = 437.569
q1a1/no_index = 0

q1a2/min_fs = 388
q1a2/min_ss = 185
q1a2/max_fs = 581
q1a2/max_ss = 369
q1a2/fs = -0.999984x -0.005516y
q1a2/ss = +0.005516x -0.999984y
q1a2/corner_x = 41.8532
q1a2/corner_y = 229.946
q1a2/no_index = 0

q1a3/min_fs = 582
q1a3/min_ss = 185
q1a3/max_fs = 775
q1a3/max_ss = 369
q1a3/fs = -0.999984x -0.005516y
q1a3/ss = +0.005516x -0.999984y
q1a3/corner_x = -155.144
q1a3/corner_y = 228.858
q1a3/no_index = 0

q1a4/min_fs = 388
q1a4/min_ss = 370
q1a4/max_fs = 581
q1a4/max_ss = 554
q1a4/fs = +0.003839x -0.999993y
q1a4/ss = +0.999993x +0.003839y
q1a4/corner_x = -355.118
q1a4/corner_y = 856.777
q1a4/no_index = 0

q1a5/min_fs = 582
q1a5/min_ss = 370
q1a5/max_fs = 775
q1a5/max_ss = 554
q1a5/fs = +0.003839x -0.999993y
q1a5/ss = +0.999993x +0.003839y
q1a5/corner_x = -354.361
q1a5/corner_y = 659.779
q1a5/no_index = 0

q1a6/min_fs = 388
q1a6/min_ss = 555
q1a6/max_fs = 581
q1a6/max_ss = 739
q1a6/fs = +0.004350x -0.999991y
q1a6/ss = +0.999991x +0.004350y
q1a6/corner_x = -147.666
q1a6/corner_y = 857.808
q1a6/no_index = 0

q1a7/min_fs = 582
q1a7/min_ss = 555
q1a7/max_fs = 775
q1a7/max_ss = 739
q1a7/fs = +0.004350x -0.999991y
q1a7/ss = +0.999991x +0.004350y
q1a7/corner_x = -146.81
q1a7/corner_y = 660.81
q1a7/no_index = 0

q1a8/min_fs = 388
q1a8/min_ss = 740
q1a8/max_fs = 581
q1a8/max_ss = 924
q1a8/fs = +0.999993x +0.003595y
q1a8/ss = -0.003595x +0.999993y
q1a8/corner_x = -781.62
q1a8/corner_y = 463.805
q1a8/no_index = 0

q1a9/min_fs = 582
q1a9/min_ss = 740
q1a9/max_fs = 775
q1a9/max_ss = 924
q1a9/fs = +0.999993x +0.003595y
q1a9/ss = -0.003595x +0.999993y
q1a9/corner_x = -584.622
q1a9/corner_y = 464.512
q1a9/no_index = 0

q1a10/min_fs = 388
q1a10/min_ss = 925
q1a10/max_fs = 581
q1a10/max_ss = 1109
q1a10/fs = +0.999977x +0.006752y
q1a10/ss = -0.006752x +0.999977y
q1a10/corner_x = -782.559
q1a10/corner_y = 673.949
q1a10/no_index = 0

q1a11/min_fs = 582
q1a11/min_ss = 925
q1a11/max_fs = 775
q1a11/max_ss = 1109
q1a11/fs = +0.999977x +0.006752y
q1a11/ss = -0.006752x +0.999977y
q1a11/corner_x = -585.563
q1a11/corner_y = 675.278
q1a11/no_index = 0

q1a12/min_fs = 388
q1a12/min_ss = 1110
q1a12/max_fs = 581
q1a12/max_ss = 1294
q1a12/fs = +0.004489x -0.999990y
q1a12/ss = +0.999990x +0.004489y
q1a12/corner_x = -761.306
q1a12/corner_y = 432.29
q1a12/no_index = 0

q1a13/min_fs = 582
q1a13/min_ss = 1110
q1a13/max_fs = 775
q1a13/max_ss = 1294
q1a13/fs = +0.004489x -0.999990y
q1a13/ss = +0.999990x +0.004489y
q1a13/corner_x = -760.423
q1a13/corner_y = 235.291
q1a13/no_index = 0

q1a14/min_fs = 388
q1a14/min_ss = 1295
q1a14/max_fs = 581
q1a14/max_ss = 1479
q1a14/fs = +0.003218x -0.999994y
q1a14/ss = +0.999994x +0.003218y
q1a14/corner_x = -546.114
q1a14/corner_y = 433.719
q1a14/no_index = 0

q1a15/min_fs = 582
q1a15/min_ss = 1295
q1a15/max_fs = 775
q1a15/max_ss = 1479
q1a15/fs = +0.003218x -0.999994y
q1a15/ss = +0.999994x +0.003218y
q1a15/corner_x = -545.478
q1a15/corner_y = 236.719
q1a15/no_index = 0

q2a0/min_fs = 776
q2a0/min_ss = 0
q2a0/max_fs = 969
q2a0/max_ss = 184
q2a0/fs = +0.008577x -0.999963y
q2a0/ss = +0.999963x +0.008577y
q2a0/corner_x = -430.54
q2a0/corner_y = 26.2725
q2a0/no_index = 0

q2a1/min_fs = 970
q2a1/min_ss = 0
q2a1/max_fs = 1163
q2a1/max_ss = 184
q2a1/fs = +0.008577x -0.999963y
q2a1/ss = +0.999963x +0.008577y
q2a1/corner_x = -428.851
q2a1/corner_y = -170.72
q2a1/no_index = 0

q2a2/min_fs = 776
q2a2/min_ss = 185
q2a2/max_fs = 969
q2a2/max_ss = 369
q2a2/fs = +0.005987x -0.999981y
q2a2/ss = +0.999981x +0.005987y
q2a2/corner_x = -223.206
q2a2/corner_y = 27.6567
q2a2/no_index = 0

q2a3/min_fs = 970
q2a3/min_ss = 185
q2a3/max_fs = 1163
q2a3/max_ss = 369
q2a3/fs = +0.005987x -0.999981y
q2a3/ss = +0.999981x +0.005987y
q2a3/corner_x = -222.026
q2a3/corner_y = -169.34
q2a3/no_index = 0

q2a4/min_fs = 776
q2a4/min_ss = 370
q2a4/max_fs = 969
q2a4/max_ss = 554
q2a4/fs = +0.999980x +0.006459y
q2a4/ss = -0.006459x +0.999980y
q2a4/corner_x = -853.405
q2a4/corner_y = -366.489
q2a4/no_index = 0

q2a5/min_fs = 970
q2a5/min_ss = 370
q2a5/max_fs = 1163
q2a5/max_ss = 554
q2a5/fs = +0.999980x +0.006459y
q2a5/ss = -0.006459x +0.999980y
q2a5/corner_x = -656.41
q2a5/corner_y = -365.218
q2a5/no_index = 0

q2a6/min_fs = 776
q2a6/min_ss = 555
q2a6/max_fs = 969
q2a6/max_ss = 739
q2a6/fs = +0.999976x +0.007031y
q2a6/ss = -0.007031x +0.999976y
q2a6/corner_x = -854.029
q2a6/corner_y = -162.264
q2a6/no_index = 0

q2a7/min_fs = 970
q2a7/min_ss = 555
q2a7/max_fs = 1163
q2a7/max_ss = 739
q2a7/fs = +0.999976x +0.007031y
q2a7/ss = -0.007031x +0.999976y
q2a7/corner_x = -657.033
q2a7/corner_y = -160.88
q2a7/no_index = 0

q2a8/min_fs = 776
q2a8/min_ss = 740
q2a8/max_fs = 969
q2a8/max_ss = 924
q2a8/fs = -0.006619x +0.999977y
q2a8/ss = -0.999977x -0.006619y
q2a8/corner_x = -458.301
q2a8/corner_y = -794.931
q2a8/no_index = 0

q2a9/min_fs = 970
q2a9/min_ss = 740
q2a9/max_fs = 1163
q2a9/max_ss = 924
q2a9/fs = -0.006619x +0.999977y
q2a9/ss = -0.999977x -0.006619y
q2a9/corner_x = -459.605
q2a9/corner_y = -597.935
q2a9/no_index = 0

q2a10/min_fs = 776
q2a10/min_ss = 925
q2a10/max_fs = 969
q2a10/max_ss = 1109
q2a10/fs = -0.009776x +0.999951y
q2a10/ss = -0.999951x -0.009776y
q2a10/corner_x = -668.447
q2a10/corner_y = -795.385
q2a10/no_index = 0

q2a11/min_fs = 970
q2a11/min_ss = 925
q2a11/max_fs = 1163
q2a11/max_ss = 1109
q2a11/fs = -0.009776x +0.999951y
q2a11/ss = -0.999951x -0.009776y
q2a11/corner_x = -670.372
q2a11/corner_y = -598.395
q2a11/no_index = 0

q2a12/min_fs = 776
q2a12/min_ss = 1110
q2a12/max_fs = 969
q2a12/max_ss = 1294
q2a12/fs = +0.999975x +0.007077y
q2a12/ss = -0.007077x +0.999975y
q2a12/corner_x = -422.972
q2a12/corner_y = -768.929
q2a12/no_index = 0

q2a13/min_fs = 970
q2a13/min_ss = 1110
q2a13/max_fs = 1163
q2a13/max_ss = 1294
q2a13/fs = +0.999975x +0.007077y
q2a13/ss = -0.007077x +0.999975y
q2a13/corner_x = -225.979
q2a13/corner_y = -767.534
q2a13/no_index = 0

q2a14/min_fs = 776
q2a14/min_ss = 1295
q2a14/max_fs = 969
q2a14/max_ss = 1479
q2a14/fs = +0.999977x +0.006867y
q2a14/ss = -0.006867x +0.999977y
q2a14/corner_x = -424.608
q2a14/corner_y = -561.775
q2a14/no_index = 0

q2a15/min_fs = 970
q2a15/min_ss = 1295
q2a15/max_fs = 1163
q2a15/max_ss = 1479
q2a15/fs = +0.999977x +0.006867y
q2a15/ss = -0.006867x +0.999977y
q2a15/corner_x = -227.611
q2a15/corner_y = -560.423
q2a15/no_index = 0

q3a0/min_fs = 1164
q3a0/min_ss = 0
q3a0/max_fs = 1357
q3a0/max_ss = 184
q3a0/fs = +0.999986x +0.005156y
q3a0/ss = -0.005156x +0.999986y
q3a0/corner_x = -21.762
q3a0/corner_y = -443.157
q3a0/no_index = 0

q3a1/min_fs = 1358
q3a1/min_ss = 0
q3a1/max_fs = 1551
q3a1/max_ss = 184
q3a1/fs = +0.999986x +0.005156y
q3a1/ss = -0.005156x +0.999986y
q3a1/corner_x = 175.235
q3a1/corner_y = -442.141
q3a1/no_index = 0

q3a2/min_fs = 1164
q3a2/min_ss = 185
q3a2/max_fs = 1357
q3a2/max_ss = 369
q3a2/fs = +0.999976x +0.007050y
q3a2/ss = -0.007050x +0.999976y
q3a2/corner_x = -22.9442
q3a2/corner_y = -238.237
q3a2/no_index = 0

q3a3/min_fs = 1358
q3a3/min_ss = 185
q3a3/max_fs = 1551
q3a3/max_ss = 369
q3a3/fs = +0.999976x +0.007050y
q3a3/ss = -0.007050x +0.999976y
q3a3/corner_x = 174.051
q3a3/corner_y = -236.848
q3a3/no_index = 0

q3a4/min_fs = 1164
q3a4/min_ss = 370
q3a4/max_fs = 1357
q3a4/max_ss = 554
q3a4/fs = -0.000308x +1.000001y
q3a4/ss = -1.000001x -0.000308y
q3a4/corner_x = 371.138
q3a4/corner_y = -866.503
q3a4/no_index = 0

q3a5/min_fs = 1358
q3a5/min_ss = 370
q3a5/max_fs = 1551
q3a5/max_ss = 554
q3a5/fs = -0.000308x +1.000001y
q3a5/ss = -1.000001x -0.000308y
q3a5/corner_x = 371.076
q3a5/corner_y = -669.505
q3a5/no_index = 0

q3a6/min_fs = 1164
q3a6/min_ss = 555
q3a6/max_fs = 1357
q3a6/max_ss = 739
q3a6/fs = -0.006246x +0.999981y
q3a6/ss = -0.999981x -0.006246y
q3a6/corner_x = 169.579
q3a6/corner_y = -866.748
q3a6/no_index = 0

q3a7/min_fs = 1358
q3a7/min_ss = 555
q3a7/max_fs = 1551
q3a7/max_ss = 739
q3a7/fs = -0.006246x +0.999981y
q3a7/ss = -0.999981x -0.006246y
q3a7/corner_x = 168.348
q3a7/corner_y = -669.753
q3a7/no_index = 0

q3a8/min_fs = 1164
q3a8/min_ss = 740
q3a8/max_fs = 1357
q3a8/max_ss = 924
q3a8/fs = -0.999997x +0.002389y
q3a8/ss = -0.002389x -0.999997y
q3a8/corner_x = 795.262
q3a8/corner_y = -472.626
q3a8/no_index = 0

q3a9/min_fs = 1358
q3a9/min_ss = 740
q3a9/max_fs = 1551
q3a9/max_ss = 924
q3a9/fs = -0.999997x +0.002389y
q3a9/ss = -0.002389x -0.999997y
q3a9/corner_x = 598.263
q3a9/corner_y = -472.157
q3a9/no_index = 0

q3a10/min_fs = 1164
q3a10/min_ss = 925
q3a10/max_fs = 1357
q3a10/max_ss = 1109
q3a10/fs = -0.999979x -0.006469y
q3a10/ss = +0.006469x -0.999979y
q3a10/corner_x = 800.574
q3a10/corner_y = -682.811
q3a10/no_index = 0

q3a11/min_fs = 1358
q3a11/min_ss = 925
q3a11/max_fs = 1551
q3a11/max_ss = 1109
q3a11/fs = -0.999979x -0.006469y
q3a11/ss = +0.006469x -0.999979y
q3a11/corner_x = 603.577
q3a11/corner_y = -684.087
q3a11/no_index = 0

q3a12/min_fs = 1164
q3a12/min_ss = 1110
q3a12/max_fs = 1357
q3a12/max_ss = 1294
q3a12/fs = -0.004206x +0.999991y
q3a12/ss = -0.999991x -0.004206y
q3a12/corner_x = 779.821
q3a12/corner_y = -441.11
q3a12/no_index = 0

q3a13/min_fs = 1358
q3a13/min_ss = 1110
q3a13/max_fs = 1551
q3a13/max_ss = 1294
q3a13/fs = -0.004206x +0.999991y
q3a13/ss = -0.999991x -0.004206y
q3a13/corner_x = 778.99
q3a13/corner_y = -244.113
q3a13/no_index = 0

q3a14/min_fs = 1164
q3a14/min_ss = 1295
q3a14/max_fs = 1357
q3a14/max_ss = 1479
q3a14/fs = -0.004897x +0.999988y
q3a14/ss = -0.999988x -0.004897y
q3a14/corner_x = 565.445
q3a14/corner_y = -439.944
q3a14/no_index = 0

q3a15/min_fs = 1358
q3a15/min_ss = 1295
q3a15/max_fs = 1551
q3a15/max_ss = 1479
q3a15/fs = -0.004897x +0.999988y
q3a15/ss = -0.999988x -0.004897y
q3a15/corner_x = 564.481
q3a15/corner_y = -242.944
q3a15/no_index = 0



q0a0/coffset = 0.567409
q0a1/coffset = 0.567409
q0a2/coffset = 0.567409
q0a3/coffset = 0.567409
q0a4/coffset = 0.567409
q0a5/coffset = 0.567409
q0a6/coffset = 0.567409
q0a7/coffset = 0.567409
q0a8/coffset = 0.567409
q0a9/coffset = 0.567409
q0a10/coffset = 0.567409
q0a11/coffset = 0.567409
q0a12/coffset = 0.567409
q0a13/coffset = 0.567409
q0a14/coffset = 0.567409
q0a15/coffset = 0.567409
q1a0/coffset = 0.567409
q1a1/coffset = 0.567409
q1a2/coffset = 0.567409
q1a3/coffset = 0.567409
q1a4/coffset = 0.567409
q1a5/coffset = 0.567409
q1a6/coffset = 0.567409
q1a7/coffset = 0.567409
q1a8/coffset = 0.567409
q1a9/coffset = 0.567409
q1a10/coffset = 0.567409
q1a11/coffset = 0.567409
q1a12/coffset = 0.567409
q1a13/coffset = 0.567409
q1a14/coffset = 0.567409
q1a15/coffset = 0.567409
q2a0/coffset = 0.567409
q2a1/coffset = 0.567409
q2a2/coffset = 0.567409
q2a3/coffset = 0.567409
q2a4/coffset = 0.567409
q2a5/coffset = 0.567409
q2a6/coffset = 0.567409
q2a7/coffset = 0.567409
q2a8/coffset = 0.567409
q2a9/coffset = 0.567409
q2a10/coffset = 0.567409
q2a11/coffset = 0.567409
q2a12/coffset = 0.567409
q2a13/coffset = 0.567409
q2a14/coffset = 0.567409
q2a15/coffset = 0.567409
q3a0/coffset = 0.567409
q3a1/coffset = 0.567409
q3a2/coffset = 0.567409
q3a3/coffset = 0.567409
q3a4/coffset = 0.567409
q3a5/coffset = 0.567409
q3a6/coffset = 0.567409
q3a7/coffset = 0.567409
q3a8/coffset = 0.567409
q3a9/coffset = 0.567409
q3a10/coffset = 0.567409
q3a11/coffset = 0.567409
q3a12/coffset = 0.567409
q3a13/coffset = 0.567409
q3a14/coffset = 0.567409
q3a15/coffset = 0.567409
