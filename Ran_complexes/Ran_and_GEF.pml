load "/Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_complexes/1i2m.pdb"
bg_color white
create Ran, chain A
create GEF, chain B
delete 1i2m
set_view (\
     0.509908617,   -0.118909307,    0.851966500,\
    -0.855097055,    0.037891444,    0.517073452,\
    -0.093766958,   -0.992177963,   -0.082355812,\
     0.000000000,    0.000000000, -178.053451538,\
     7.465000629,  -16.553627014,    5.558338642,\
   140.378692627,  215.728210449,  -20.000000000 )
alter open_conf, chain = "A"
alter all, segi = ""
color density, Ran  
color deepteal, GEF
select P_loop, Ran and resi 17-23
color red, P_loop
select switch1, Ran and resi 39-45
select switch2, Ran and resi 67-75
color limon, switch1
color yelloworange, switch2
select R108, Ran and resi 106
select K101, Ran and resi 99
select R78, Ran and resi 76
select D79, Ran and resi 77
select R112, Ran and resi 110
color gray70, not polymer
as cartoon
show sticks, not polymer
show sticks, R108 or K101 or R78 or D79 or R112
color hotpink, R108 or K101 or R78 or D79 or R112
select mutations, Ran and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178
