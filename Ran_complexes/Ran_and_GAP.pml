load "/Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_complexes/1k5d.pdb"
bg_color white
create Ran, chain A
create YRB1, chain B
create GAP, chain C
delete 1k5d
set_view (\
     0.261730522,    0.360086948,    0.895447671,\
    -0.195272952,   -0.888845265,    0.414507449,\
     0.945170820,   -0.283344537,   -0.162326813,\
     0.000000000,    0.000000000, -210.501907349,\
     2.396888733,   -2.288713455,  -33.023170471,\
   165.961303711,  255.042510986,  -20.000000000 )
alter open_conf, chain = "A"
alter all, segi = ""
color density, Ran  
color deeppurple, GAP
color yelloworange, YRB1
select P_loop, Ran and resi 17-23
color red, P_loop
select switch1, Ran and resi 39-45
select switch2, Ran and resi 67-75
color limon, switch1
color yelloworange, switch2
select T34, Ran and resi 32
select Q147, Ran and resi 145
select F58, Ran and resi 56
select D79, Ran and resi 77
select Y157, Ran and resi 155
select K132, Ran and resi 130
color gray70, not polymer
as cartoon
show sticks, not polymer
show spheres, T34 or Q147 or F58 or Y157 or K132 or D79
color hotpink, T34 or Q147 or F58 or Y157 or K132 or D79
select mutations, Ran and resi 32+56+76+77+78+82+99+100+103+106+110+113+127+130+135+137+139+141+145+146+152+155+167+178
