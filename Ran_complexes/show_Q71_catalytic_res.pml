load "/Users/tperica/Documents/Gsp1_bioinformatics/Ran_structures/Ran_complexes/1k5d.pdb"
bg_color white
create Ran, chain A
create YRB1, chain B
create GAP, chain C
select tail, Ran and resi 179-213
delete 1k5d or YRB1 or GAP or tail
select switch2, resi 67-75
select switch1, resi 39-45
select Q71, resi 69
color density, Ran  
show sticks, not polymer
zoom vis
