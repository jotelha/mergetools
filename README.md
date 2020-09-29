# VMD MergeTools package. Version 0.1


## Overview

MergeTools is a plugin for [VMD](http://www.ks.uiuc.edu/Research/vmd/)
that helps to merge two (solvated) systems.

## Dependencies

- required:
  - topotools 1.8
- optional:
  - pbctools 3.0:
    - tcllib 1.20

## Interface

After importing via `package require mergetools`, all commands are accessed via
the command `mrg`. In analogy to `topotools` and other VMD commands, `mrg`
always accepts the `-molid` and `-sel` options. Many subcommands follow the
pattern

    mrg -sel $base $subcommand (and possibly $subsubcommand) $ext

where `$base` and `$ext` are atom selections.

Valid `$subommand`s and `$subsubcommand`s are

* `set` sets MergeTools-wide properties
  * `autowrap`, defaults to `atom`, valid choices are `atom`, `res`, `resid`,
                `residue`, `seg`, `segid`, `chain`, `fragment`, `off`, (de-)
                activates automatic wrapping in `merge` command.
  * `autojoin`, defaults to `residue`,
                valid choices are `res`, `resid`, `residue`, `seg`, `segid`,
                `chain`, `fragment`, `bonded`, `connected`, `off`, (de-)
                activates automatic wrapping in `merge` command. First
                wrapping and then joining can fix "dirty" systems where
                molecules might be split across periodic boundaries without
                proper track-keeping of image flags.
  * `bondlist`, defaults to `on`, valid choces are `on` and `off`. Uses
                reliable, but slow algorithm for joining split compounds.
                See the `pbctools` (https://github.com/frobnitzem/pbctools)
                documentation for details.
  * `compound`, defaults to `residue`, modify to operate on different hierarchy
                level of compounds.
  * `compname`, defaults to `resname`, modify to operate on different hierarchy
                level of compounds.
  * `dispensable`, `SOL` by default.
  * `forbidden`, not `dispensable` by default.
  * `immobile`, None by default.
  * `mobile`, not `dispensable | immobile` by default.
  * `overlap`, defaults to 2.0 (Angstrom). Distance tolerance below which two
               atoms are identified as overlapping.
  * `maxit`, defaults to 1000. Maximum number of iterations for random
             movements of a single component within the `move` subcommand.
* `get` is used to query above properties and in addition
  * `base` and
  * `ext`, the atom selections that arise in a `merge` command call, as well as
  * `compnames`, list of unique component names in molecule or selection.
* `merge`, merges two selections from systems of different `molid`, i.e.
           `mrg -sel $base merge $ext`, where a new system with distinct
           `molid` and same cell as `base` arises. Returns new `molid`-
* `overlap`, identifies several overlaps, usually between selections `base` and
            `ext` in the same system, i.e.
  * `mobile`, all mobile compounds in `ext` that overlap with forbidden
              compounds in `base` and vice versa.
  * `removable`, all dispensable compounds in `ext` that overlap with
                 compounds in `base` and remaining dispensable compounds in
                 `base` that overlap with not yet selected compounds in `ext`.
  * `dispensable`, all dispensable compounds in `ext` that overlap with any
                   compound in `base` and vice versa.
  * `forbidden`, all compounds in `ext` that overlap with forbidden compounds
                 in `base` and vice versa.
  * `unidirectional`, all compounds in `ext` that overlap with any compound in `base`.
  * `bidirectional`, all compounds in `ext` that overlap with any compound in
                     `base` and vice versa.
  * `internal`, only accepts `base`. Identifies internal overlap.
                Experimental, untested and slow.
* `move`, move every compound in `ext` around randomly until it has no forbidden
          overlap with `base` anymore. Positions are randomly generated within
          the bounding box of `base`. Alternatively, a bounding box for random
          coordinates can be specified as second positional argument, i.e.
          `mrg -sel $base move $ext [[$xlo $ylo $zlo][$xhi $yhi $zhi]]`.
* `remove`, creates new system from `base` under removal of `ext`.
            Returns new `molid`-
* `report`
  * `compounds`, prints summary of categorized compounds in `base`.
  * `overlap`, prints summary of categorized overlap of `ext` with `base`.
* `help`, print short usage overview.

## Sample usage

When merging solvated systems, i.e. a gold nanoparticle above a gold substrate
in aqueous electrolyte solution, overlapping compounds must be dealt with. When
handling overlap in merged system, MergeTools distinguishes four categories of
compounds: immobile, mobile, dispensable and forbidden. In our example,
nanoparticle and substrate will constitute the immobile compounds of the system,
as we want to align them as we desire without MergeTools touching them. The
solvent is dispensable, meaning that MergeTools is allowed to remove overlapping
water molecules. Electrolyte ions must not be removed arbitrarily in order to
conserve the system's total charge, usually zero. However, they might be allowed
to be moved around a bit without disturbing a system's (equilibrated) state too
much. Thus they are marked as mobile and may be assigned random positions in
case of overlap. The forbidden selection is the negation of the dispensable
selection by default, meaning that any overlap conflict involving forbidden
compounds must be resolved by moving mobile or removing dispensable compounds.

```tcl
package require topotools
package require pbctools
package require mergetools

# read files
set substrate_system_id [mol new substrate.gro waitfor all]
set particle_system_id [mol new particle.gro waitfor all]
mol top $substrate_system_id

# substrate system selections
set substrate_system [atomselect $substrate_system_id all]

# particle system selections
set particle_system [atomselect $particle_system_id all]

# merge command allows to automatically wrap and join via pbctools if desired
mrg set autowrap residue autojoin residue
mrg set compound residue compname resname
mrg set dispensable SOL
mrg set immobile AUM
# remaining compounds are treated as 'mobile' if not set explicitly
mrg set mobile NA CL

# list categorized compunds in both systems
mrg -sel $substrate_system report compounds
mrg -sel $particle_system report compounds

# mrg merge leaves the two original systems untouched and uses topotools under
# the hood to create a new merged system whose molid is returned:
set merged_id [mrg -sel $substrate_system merge $particle_system]
# First and second atom selection must exist within systems of different molid.
# The first selection as an argument to '-sel' is treated as the 'base' system
# and the merged system will inherit its cell. The positional argument to
# 'merge' is 'ext' (extension) system. The distinction between those will play
# a subtle role for most subsequent commands.

# Selections representing the original two systems within the newly created
# merged system are obtained via
set base [mrg get base]
set ext [mrg get ext]
# where 'base' will represent 'substrate' and 'ext' particle' in this example.
# Subsequently, MergeTools makes use of those two disjoint selections to
# identify overlap efficiently. The core assumption here is that no internal
# overlap exists within the two constituent systems.

set merged [atomselect $merged_id all]
mrg -sel $merged report compounds

# The 'overlap mobile' subcommand yields all 'mobile' compounds that overlap
# with 'forbidden' compounds.
set overlap_to_move [mrg -sel $base overlap mobile $ext]
# This overlap conflict can be alleviated by moving all overlapping 'mobile'
# compounds to allowed locations, i.e. locations where they would overlap
# with 'dispensable' compounds at most.
mrg -sel $merged move $overlap_to_move

# Remaining overlap conflicts are resolved by identifying a certain set of
# 'dispensable' compounds whose removal will render the remaining system
# overlap-free (if possible).
set overlap_to_remove [mrg -sel $base overlap removable $ext]
# The removal strategy first selects all 'dispensable' compounds in 'ext' that
# overlap with anything in 'base', then selects 'dispensable' compounds in
# 'base' that overlap with anything in 'ext' that has not been selected.
# Just as 'merge' will 'remove' leave the input system attached and elevate
# the system subtracted by the identified dispensable overlap to a new molid.
set nonoverlap_id [mrg -sel $merged remove $overlap_to_remove]
mrg -molid $nonoverlap_id -sel all report compound

# visualization
pbc box -on

mol off $substrate_system_id
mol off $particle_system_id
mol off $merged_id
mol on $nonoverlap_id

# merged systen representations
mol delrep 0 $merged_id
set nonsolvent_rep 0
set base_rep 1
set ext_rep  2
set removed_overlap_rep  3

for {set i 0} {$i < 4} {incr i} {
    mol addrep $merged_id
}

mol modselect $nonsolvent_rep $merged_id "not resname SOL"
mol modselect $base_rep $merged_id [$base text]
mol modselect $ext_rep $merged_id [$ext text]
mol modselect $removed_overlap_rep $merged_id [$overlap_to_remove text]

mol showrep $merged_id $base_rep off
mol showrep $merged_id $ext_rep off
mol showrep $merged_id $nonsolvent_rep off
mol showrep $merged_id $removed_overlap_rep on

# overlap-free system representations
mol modselect 0 $nonoverlap_id "not resname SOL"
```

## Installation

MergeTools is written entirely in the Tcl scripting language for use as a plugin
with VMD. You can place all files within the same directory as this `README.md`
file within your VMD Plugin directory `$PLUGINDIR` under
`$PLUGINDIR/noarch/tcl/mergetools${MRG_VERSION}`, where `${MRG_VERSION}` must
be the current version of this package, i.e. `0.1`. Alternatively, place it at
any location (still within a directory `mergetools${MRG_VERSION}`) and make
sure to point the environment variable `TCLLIBPATH` to the parent directory
of this location before launching VMD, i.e.
`export TCLLIBPATH="$HOME/path/to/my/tcl/libs $TCLLIBPATH"` if this plugin is
placed at `$HOME/path/to/my/tcl/libs/mergetools0.1`. Note that the path
separator for `TCLLIBPATH` is a single space -- other than the colon typical
for path variables on UNIX systems.

## Acknowledgements

Code structure and interface based upon

    TopoTools v1.8 by Axel Kohlmeyer & Josh Vermaas, 2019.
    (https://doi.org/10.5281/zenodo.598373)

and the examples on https://sites.google.com/site/akohlmey/software/topotools/topotools-tutorial---various-tips-tricks.

## Feedback

Please report an issues at https://github.com/jotelha/mergetools/issues.
