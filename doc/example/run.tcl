# mergetools/doc/example/run.tcl
#
# Copyright (C) 2020 IMTEK Simulation
# Author: Johannes Hoermann, johannes.hoermann@imtek.uni-freiburg.de
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This script illustrates how to load, align and merge two solvated systems
# with VMD, mergetools, topotools, and pbctools. Run with
#
#    vmd -e run.tcl
#
# note: the ;list appendix to every line suppresses any undesired output
package require topotools; list
package require pbctools; list
package require mergetools; list

# tolerance for aligning substrate at box boundary
set tol 2.0; list

set z_dist 50.0; list
set x_shift 30.0; list
set y_shift 0.0; list

set solvent_resname "SOL"; list
set substrate_resname "AUM"; list
set surfactant_resname "SDS"; list
set counterion_resname "NA"; list

# read files
set substrate_system_id [mol new dat/substrate.gro waitfor all autobonds off]; list
set particle_system_id [mol new dat/particle.gro waitfor all autobonds off]; list
mol top $substrate_system_id; list

# substrate system selections
set substrate_system [atomselect $substrate_system_id all]; list
set substrate [atomselect $substrate_system_id "resname $substrate_resname"]; list

set substrate_system_bb [measure minmax $substrate_system]; list
set substrate_bb [measure minmax $substrate]; list
set substrate_system_cell [pbc get -molid $substrate_system_id]; list

# particle system selections
set particle_system [atomselect $particle_system_id all]; list
set particle [atomselect $particle_system_id "resname $substrate_resname"]; list

set particle_system_bb [measure minmax $particle_system]; list
set particle_bb [measure minmax $particle]; list
set particle_system_cell [pbc get -molid $particle_system_id]; list

# Fix resid and residue issue:
# resid holds the residue identifier as found within the datafile
#   its maximum is limited in both pdb and gro files (100000 for gro)
#   and simply starts at zero again. 
# residue is vmd's internal resiude numbering
#   and is only assigned with autobonds inference, which cannot be relied
#   on in case of residues split across PBC.
# thus, assign consecutive resid based on resid change, not absolute resid:
set resid_offsets [list]; list
set consecutive_resid 0; list
foreach sel [list $substrate_system $particle_system] {
    lappend resid_offsets $consecutive_resid; list
    set prev_resid -1; list
    set resid_list  [list]; list
    foreach cur_resid [$sel get resid] {
        if { $cur_resid != $prev_resid } {
            set prev_resid $cur_resid; list
            incr consecutive_resid; list
        }
        lappend resid_list $consecutive_resid; list
    }
    $sel set resid $resid_list; list
}

# log info
vmdcon -info [format "%-30.30s" "  #atoms in substrate system:"] \
    [format "%12d" [$substrate_system num]]
vmdcon -info [format "%-30.30s" "  #atoms in substrate:"] \
    [format "%12d" [$substrate num]]

vmdcon -info [format "%-30.30s" "  #atoms in particle system:"] \
    [format "%12d" [$particle_system num]]
vmdcon -info [format "%-30.30s" "  #atoms in particle:"] \
    [format "%12d" [$particle num]]

vmdcon -info [format "%-30.30s" "  #resid in substrate system:"] \
    [format "%12d" [llength [lsort -unique [$substrate_system get resid]]]]
vmdcon -info [format "%-30.30s" "  #resid in particle system:"] \
    [format "%12d" [llength [lsort -unique [$particle_system get resid]]]]

vmdcon -info [format "%-30.30s" "  residue names in substrate system:"] \
    [lsort -unique [$substrate_system get resname]]
vmdcon -info [format "%-30.30s" "  residue names in particle system:"] \
    [lsort -unique [$particle_system get resname]]

vmdcon -info [format "%-30.30s" "  substrate sys cell:"] $substrate_system_cell
vmdcon -info [format "%-30.30s" "  substrate sys bounding box:"] $substrate_system_bb
vmdcon -info [format "%-30.30s" "  substrate bounding box:"] $substrate_bb

vmdcon -info [format "%-30.30s" "  particle sys cell:"] $particle_system_cell
vmdcon -info [format "%-30.30s" "  particle sys bounding box:"] $particle_system_bb
vmdcon -info [format "%-30.30s" "  particle bounding box:"] $particle_bb

# z shift (align substrate at boundary)
set z_shift [ list 0. 0. [ expr [lindex $substrate_system_bb 0 2] - [lindex $substrate_bb 0 2] + $tol ] ]; list

$substrate_system moveby $z_shift; list
pbc wrap -molid $substrate_system_id -compound resid -all -verbose; list

# reset substrate_bb after z shift
set substrate_bb [measure minmax $substrate]; list

set z_shift [ list 0. 0. [ expr [lindex $substrate_bb 1 2] - [lindex $particle_bb 0 2] + $z_dist ] ]; list
$particle_system moveby $z_shift; list

# xy shift
set substrate_com [measure center $substrate]; list
set particle_com [measure center $particle]; list

set xy_shift [list \
    [ expr [lindex $substrate_com 0] - [lindex $particle_com 0] + $x_shift ] \
    [ expr [lindex $substrate_com 1] - [lindex $particle_com 1] + $y_shift ] \
    0.0 \
]; list
$particle_system moveby $xy_shift; list

## mergetools
vmdcon -info "############################"
vmdcon -info "### configure MergeTools ###"
vmdcon -info "############################"

mrg set autowrap atom autojoin resid; list
mrg set compound resid compname resname; list
mrg set bondlist on; list
mrg set dispensable $solvent_resname; list
mrg set immobile $substrate_resname; list
mrg set mobile $surfactant_resname $counterion_resname; list

vmdcon -info ""
vmdcon -info "######################################################"
vmdcon -info "### substrate and particle system compound reports ###"
vmdcon -info "######################################################"
vmdcon -info ""
mrg -sel $substrate_system report compounds; list
mrg -sel $particle_system report compounds; list

# merge
vmdcon -info ""
vmdcon -info "###########################################"
vmdcon -info "### substrate and particle system merge ###"
vmdcon -info "###########################################"
vmdcon -info ""
set merged_id [mrg -sel $substrate_system merge $particle_system]; list
set base [mrg get base]; list
set ext [mrg get ext]; list

vmdcon -info ""
vmdcon -info "#####################################"
vmdcon -info "### merged system compound report ###"
vmdcon -info "#####################################"
vmdcon -info ""
set merged [atomselect $merged_id all]; list
mrg -sel $merged report compounds; list

# move
set overlap_to_move [mrg -sel $base overlap mobile $ext]; list
vmdcon -info ""
vmdcon -info "############################################"
vmdcon -info "### merged system movable overlap report ###"
vmdcon -info "############################################"
vmdcon -info ""
mrg -sel $overlap_to_move report compounds; list
mrg -sel $merged move $overlap_to_move; list

# remove
set overlap_to_remove [mrg -sel $base overlap removable $ext]; list
vmdcon -info ""
vmdcon -info "##############################################"
vmdcon -info "### merged system overlap to remove report ###"
vmdcon -info "##############################################"
vmdcon -info ""
mrg -sel $overlap_to_remove report compounds; list
# switch autojoin off here for performance, not safe
# mrg set autojoin off; list
set keep_id [mrg -sel $merged remove $overlap_to_remove]; list
set keep [atomselect $keep_id all]; list

mrg -molid $keep_id -sel $keep report compound; list

# write output
$keep writepdb default.pdb; list

## visualization 
pbc box -on; list

mol off $substrate_system_id; list
mol off $particle_system_id; list
mol off $merged_id; list
mol on $keep_id; list

mol delrep 0 $merged_id; list
set nonsolvent_rep 0; list
set base_rep 1; list
set ext_rep  2; list
set removed_overlap_rep  3; list

for {set i 0} {$i < 4} {incr i} {
    mol addrep $merged_id; list
}
mol modselect $nonsolvent_rep $merged_id "not resname $solvent_resname"; list
mol modselect $base_rep $merged_id [$base text]; list
mol modselect $ext_rep $merged_id [$ext text]; list
mol modselect $removed_overlap_rep $merged_id [$overlap_to_remove text]; list

mol showrep $merged_id $base_rep off; list
mol showrep $merged_id $ext_rep off; list
mol showrep $merged_id $nonsolvent_rep off; list
mol showrep $merged_id $removed_overlap_rep on; list

# overlap-free system representations
mol modselect 0 $keep_id "not resname $solvent_resname"; list