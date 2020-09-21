#!/usr/bin/tclsh
# mergetool, a VMD package to merge overlapping molecular systems in  VMD with
# the help of TopoTools and PbcTools.
#
# Copyright (c) 2018-2020
#               by Johannes Hoermann <johannes.hoermann@imtek.uni-freiburg.de>
#
# $Id: mergetool.tcl,v 0.1 2020/09/19 $
#
# Sample usage for merging to systems
#
#   vmd> package require mergetool
#   vmd> set base_id [mol new base.gro waitfor all]
#   vmd> set ext_id [mol new ext.gro waitfor all]
#   vmd> set base [atomselect $base_id all]
#   vmd> set dry_ext [atomselect $ext_id "not resname SOL"]
#   vmd> merge $base $dry_ext
#
# Double-check for section
#   Info) Identify ovelap...
#   Info) #atoms in overlapping SOD:                0
#   Info) #atoms in overlapping TIP3:            9354
#   Info) #atoms in overlapping SDS:                0
# in output to make sure only water has been removed.
#
namespace eval ::MergeTools:: {
    variable version 0.1

    package require topotools
    package require pbctools

    # distance within which molecules are regarded as overlapping
    variable overlap_distance 2.0
    variable maxit 1000

    variable autowrap "-nocompound"
    variable autojoin "residue"

    variable compound "residue"
    variable compname "resname"

    variable negative_keywords {"no" "none" "off" "false"}

    variable valid_autowrap_vals [list \
        {*}$negative_keywords \
        "atom" "res" "resid" "residue" "seg" "segid" "chain" "fragment"]

    variable valid_autojoin_vals [list \
        {*}$negative_keywords \
        "seg" "segid" "res" "resid" "residue" \
        "chain" "fragment" "bonded" "connected"]

    variable valid_compound_vals [list \
        "res" "resid" "residue" "seg" "segid" "chain" "fragment"]

    variable immobile_compounds {}
    variable mobile_compounds -1
    variable dispensable_compounds {SOL}
    variable forbidden_compounds -1

    variable base
    variable ext
}

# help/usage/error message and online documentation.
proc ::MergeTools::usage {} {
    vmdcon -info ""
    vmdcon -info "MergeTools, a VMD package to merge molecular systems with"
    vmdcon -info "configurable overlap identification, displacement and removal."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "usage:  <command> \[args...\] <flags>"
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "common flags (not implemented):"
    vmdcon -info ""
    vmdcon -info "  -molid     <num>|top    molecule id (default: 'top')"
    vmdcon -info "  -sel       <selection>  atom selection function or text (default: 'all')"
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "commands:"
    vmdcon -info ""
    vmdcon -info "  help                                          prints this message"
    vmdcon -info ""
    vmdcon -info "  info                                          display system information."
    vmdcon -info "  join <key>                                    (re-)joins residues split across boundries."
    vmdcon -info "  wrap <key>                                    wrap system into one periodic image."
    vmdcon -info ""
    vmdcon -info "key - value pairs"
    vmdcon -info ""
    vmdcon -info "  join:   residue           "
    vmdcon -info "  set:    bb                bounding box ({ { float float float } { float float float } })"
    vmdcon -info "          distance          desired distance between surface and indenter (float)"
    vmdcon -info "          interfaceInfile   input LAMMPS data file of interface (str)"
    vmdcon -info "          indenterInfile    input LAMMPS data file of indenter (str)"
    vmdcon -info "          outputPrefix      output prefix prepended to all resulting files (str)"
    vmdcon -info "  show:   nonsolvent        "
    vmdcon -info "          solvent           "
    vmdcon -info "          surfactant        "
    vmdcon -info "  use:    sds               "
    vmdcon -info "          ctab              "
    vmdcon -info "  wrap:   atom              "
    vmdcon -info "          residue           "
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "sample usage:"
    vmdcon -info ""
    vmdcon -info "  vmd> package require mergetool"
    vmdcon -info "  vmd> set base_id \[mol new base.gro waitfor all\]"
    vmdcon -info "  vmd> set ext_id \[mol new ext.gro waitfor all\]"
    vmdcon -info "  vmd> set base \[atomselect \$base_id all\]"
    vmdcon -info "  vmd> set dry_ext \[atomselect \"\$ext_id 'not resname SOL\"\]"
    vmdcon -info "  vmd> merge \$base \$dry_ext"
    vmdcon -info ""
    vmdcon -info "will use parameters for SDS, merge an interfacial system's 'interface.lammps'"
    vmdcon -info "data file with an indenter's 'indenter.pdb', remove any overlap and write"
    vmdcon -info "system.lammps, system.psf, system.pdb data files as well as a tga snapshot"
    vmdcon -info "system.tga of the resulting system."
    vmdcon -info ""
    vmdcon -info ""
    vmdcon -info "Copyright (c) 2018-2020"
    vmdcon -info "              by Johannes Hoermann <johannes.hoermann@imtek.uni-freiburg.de>"
    vmdcon -info ""
    return
}

# the main frontend command.
# this takes care of all sanity checks on arguments and
# then dispatches the subcommands to the corresponding
# subroutines.
proc ::MergeTools::mrg { args } {
    variable version

    variable negative_keywords
    variable valid_autowrap_vals
    variable valid_autojoin_vals

    set molid -1
    set seltxt all
    set localsel 1
    set selmol -1

    # variable overlap_distance

    set cmd {}
    set sel {} ; # need to initialize it here for scoping

    # process generic arguments and remove them
    # from argument list.
    set newargs {}
    for {set i 0} {$i < [llength $args]} {incr i} {
        set arg [lindex $args $i]

        if {[string match -?* $arg]} {

            set val [lindex $args [expr $i+1]]

            switch -- $arg {
                -molid {
                    if {[catch {molinfo $val get name} res]} {
                        vmdcon -err "Invalid -molid argument '$val': $res"
                        return
                    }
                    set molid $val
                    if {[string equal $molid "top"]} {
                        set molid [molinfo top]
                    }
                    incr i
                }

                -sel {
                    # check if the argument to -sel is a valid atomselect command
                    if {([info commands $val] != "") && ([string equal -length 10 $val atomselect])} {
                        set localsel 0
                        set selmol [$val molid]
                        set sel $val
                    } else {
                        set localsel 1
                        set seltxt $val
                    }
                    incr i
                }

                -- break

                default {
                    vmdcon -info "default: $arg"
                }
            }
        } else {
            lappend newargs $arg
        }
    }

    if {$molid < 0} {
        set molid $selmol
    }
    if {$molid < 0} {
        set molid [molinfo top]
    }

    set retval ""
    if {[llength $newargs] > 0} {
        set cmd [lindex $newargs 0]
        set newargs [lrange $newargs 1 end]
    } else {
        set newargs {}
        set cmd help
    }

    # check whether we have a valid command.
    set validcmd {
      "merge"
      "overlap"
      "move"
      "get"
      "set"
      "report"
      "wrap" "join"
      "help" "info" }

      # valid properties.
      set validprop {
        "autojoin"
        "autowrap"
        "base"
        "compound"
        "compnames"
        "dispensable"
        "ext"
        "forbidden"
        "immobile"
        "mobile"
    }

    if {[lsearch -exact $validcmd $cmd] < 0} {
        vmdcon -err "Unknown sub-command '$cmd'"
        usage
        return
    }

    if { ![string equal $cmd help] } {
        if {($selmol >= 0) && ($selmol != $molid)} {
            vmdcon -err "Molid from selection '$selmol' does not match -molid argument '$molid'"
            return
        }
        if {$molid < 0} {
            vmdcon -err "Cannot use 'merge $cmd' without a molecule"
            return
        }

        if {$localsel} {
            # need to create a selection
            if {[catch {atomselect $molid $seltxt} sel]} {
                vmdcon -err "Problem with atom selection using '$seltxt': $sel"
                citation_reminder
                return
            }
        }
    }

    # branch out to the various subcommands
    switch -nocase -- $cmd {
        "get" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                autojoin {
                    variable autojoin
                    if {[string equal $autojoin ""]} {
                        set retval 0
                    } else {
                        set retval $autojoin
                    }
                }

                autowrap {
                    variable autowrap
                    if {[string equal $autowrap ""]} {
                        set retval 0
                    } elseif {[string equal $autowrap "-nocompound"]} {
                        set retval "atom"
                    } else {
                        set retval [lindex [split $autowrap " "] 1]
                    }
                }

                compound {
                    variable compound
                    set retval $compound
                }

                compnames {
                    set retval [ compnames $molid $sel ]
                }

                base {
                    variable base
                    set retval $base
                }

                ext {
                    variable ext
                    set retval $ext
                }

                dispensable {
                    variable dispensable_compounds
                    set retval $dispensable_compounds
                }

                immobile {
                    variable immobile_compounds
                    set retval $immobile_compounds
                }

                mobile {
                    variable mobile_compounds
                    if { $mobile_compounds == -1 } {
                        # per default, use negative selection of immobile and dispensable
                        foreach c [ mrg get compnames ] { dict set buf $c 1 }
                        foreach c [ mrg get immobile ] { dict unset buf $c }
                        foreach c [ mrg get dispensable ] { dict unset buf $c }
                        set retval [ dict keys $buf ]
                    } else {
                        # mobile compounds set explicitly
                        set retval $mobile_compounds
                    }
                }

                forbidden {
                    variable forbidden_compounds
                    if { $forbidden_compounds == -1 } {
                        # per default, just use negation of dispensable
                        foreach c [ mrg get compnames ] { dict set buf $c 1 }
                        foreach c [ mrg get dispensable ] { dict unset buf $c }
                        set retval [ dict keys $buf ]
                    } else {
                        # forbidden compounds set explicitly
                        set retval $forbidden_compounds
                    }
                }

                default {
                    vmdcon -err "Unknown property: $key. Allowed values: $validprop"
                }
            }
        }
        "set" {
            while {[llength $newargs] > 1} {
                set key [lindex $newargs 0]
                set newargs [lrange $newargs 1 end]
                switch -nocase -- $key {
                    autojoin {
                        set val [lindex $newargs 0]
                        if {$val ni $valid_autojoin_vals} {
                        		vmdcon -err "'merge set $key' allows choice from '$valid_autojoin_vals'"
                            return
                        }
                        if {$val in $negative_keywords} {
                            variable autojoin ""
                        } else {
                            variable autojoin $val
                        }
                        set newargs [lrange $newargs 1 end]
                    }

                    autowrap {
                        set val [lindex $newargs 0]
                        if {$val ni $valid_autowrap_vals} {
                            vmdcon -err "'merge set $key' allows choice from '$valid_autowrap_vals'"
                            return
                        }
                        if {$val in $negative_keywords} {
                            variable autowrap ""
                        } elseif {[string equal $val "atom"]} {
                            variable autowrap "-nocompound"
                        } else {
                            variable autowrap "-compound $val"
                        }
                        set newargs [lrange $newargs 1 end]
                    }

                    compound {
                        set val [lindex $newargs 0]
                        if {$val ni $valid_compound_vals} {
                        		vmdcon -err "'merge set $key' allows choice from '$valid_autojoin_vals'"
                            return
                        }
                        variable compound $val
                        set newargs [lrange $newargs 1 end]
                    }

                    base {
                        variable base $sel
                    }

                    ext {
                        variable ext $sel
                    }

                    mobile {
                        if {[llength $newargs] == 0} {
                            # explicitly set empty set
                            variable mobile_compounds {}
                        } elseif {[lindex $newargs 0] in $negative_keywords} {
                            # use negative keywords to (re-)activate default behavior
                            variable mobile_compounds -1
                            set newargs [lrange $newargs 1 end]
                        } else {
                            # otherwise just jet explicitly, but check for valid entries
                            set valid_vals [ compnames $molid $sel ]
                            foreach val $newargs {
                                if {$val ni $valid_vals} {
                                    variable compound
                                    variable compname
                                    vmdcon -err "'No $compound with $compname '$val' in '[$sel text]'. Allowed choices: '$valid_vals'"
                                    return
                                }
                            }
                            variable mobile_compounds $newargs
                            set newargs {}
                        }
                    }

                    immobile {
                        if {[llength $newargs] == 0} {
                            # explicitly set empty set
                            variable immobile_compounds {}
                        } elseif {[lindex $newargs 0] in $negative_keywords} {
                            # use negative keywords in the same manner
                            variable immobile_compounds {}
                            set newargs [lrange $newargs 1 end]
                        } else {
                            # otherwise just set explicitly, but check for valid entries
                            set valid_vals [ compnames $molid $sel ]
                            foreach val $newargs {
                                if {$val ni $valid_vals} {
                                    variable compound
                                    variable compname
                                    vmdcon -err "'No $compound with $compname '$val' in '[$sel text]'. Allowed choices: '$valid_vals'"
                                    return
                                }
                            }
                            variable immobile_compounds $newargs
                            set newargs {}
                        }
                    }

                    dispensable {
                        if {[llength $newargs] == 0} {
                            # explicitly set empty set
                            variable dispensable_compounds {}
                        } elseif {[lindex $newargs 0] in $negative_keywords} {
                            # use negative keywords in the same manner
                            variable dispensable_compounds {}
                            set newargs [lrange $newargs 1 end]
                        } else {
                            # otherwise just set explicitly, but check for valid entries
                            set valid_vals [ compounds $molid $sel ]
                            foreach val $newargs {
                                if {$val ni $valid_vals} {
                                    variable compound
                                    variable compname
                                    vmdcon -err "'No $compound with $compname '$val' in '[$sel text]'. Allowed choices: '$valid_vals'"
                                    return
                                }
                            }
                            variable dispensable_compounds $newargs
                            set newargs {}
                        }
                    }

                    forbidden {
                        if {[llength $newargs] == 0} {
                            # explicitly set empty set
                            variable forbidden_compounds {}
                        } elseif {[lindex $newargs 0] in $negative_keywords} {
                            # use negative keywords to (re-)activate default behavior
                            variable forbidden_compounds -1
                            set newargs [lrange $newargs 1 end]
                        } else {
                            # otherwise just jet explicitly, but check for valid entries
                            set valid_vals [ compounds $molid $sel ]
                            foreach val $newargs {
                                if {$val ni $valid_vals} {
                                    variable compound
                                    variable compname
                                    vmdcon -err "'No $compound with $compname '$val' in '[$sel text]'. Allowed choices: '$valid_vals'"
                                    return
                                }
                            }
                            variable forbidden_compounds $newargs
                            set newargs {}
                        }
                    }

                    default {
                        vmdcon -warn "Unknown property: $key. Allowed values: $validprop"
                    }
                }
            }
            set retval 0
        }
        "merge" {
            set retval [merge $molid $sel [lindex $newargs 0]]
        }
        "overlap" {
            set retval [overlap $molid $sel [lindex $newargs 0]]
        }
        "move" {
            set retval [move $molid $sel [lindex $newargs 0]]
        }
        "report" {
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                compound -
                compounds {
                    report_compounds $molid $sel
                }

                overlap {
                    report $molid $sel [lindex $newargs 0]
                }
            }
            set retval 0
        }

        "help" {
            usage
            set retval 0
        }

        default {
            vmdcon -err "Unknown sub-command: $cmd"
            set retval 1
        }
    }
    return $retval
}

proc ::MergeTools::compnames { molid sel } {
    variable compname
    set sel [ validate_atomselect $molid $sel ]
    set retval [ lsort -unique [ $sel get $compname ] ]
}

proc ::MergeTools::validate_atomselect { molid val } {
    set seltxt all
    set localsel 1
    set selmol -1

    # check if the sel arguments are valid atomselect commands or just text
    if {([info commands $val] != "") && ([string equal -length 10 $val atomselect])} {
        set localsel 0
        set selmol [$val molid]
        set sel $val
    } else {
        set localsel 1
        set seltxt $val
    }
    if {($selmol >= 0) && ($selmol != $molid)} {
        vmdcon -err "Molid from selection '$selmol' does not match molid argument '$molid'"
        return
    }
    if {$molid < 0} {
        vmdcon -err "Cannot use 'merge $cmd' without a molecule"
        return
    }
    if {$localsel} {
        # need to create a selection
        if {[catch {atomselect $molid $seltxt} sel]} {
            vmdcon -err "Problem with atom selection using '$seltxt': $sel"
            return
        }
    }
    $sel global
    return $sel
}


proc ::MergeTools::merge { molid base_sel ext_sel } {
    variable autowrap
    variable autojoin
    variable compound
    variable compname

    # TODO: separate validation
    # set base_sel [validate_atomselect $molid $base_sel]

    set merged_id [::TopoTools::selections2mol "$base_sel $ext_sel"]
    set base_id [$base_sel molid]
    set ext_id [$ext_sel molid]

    vmdcon -info "  Merged extension 'atomselect $ext_id \"[$ext_sel text]\"' into base 'atomselect \"$base_id [$base_sel text]\"' under molid $merged_id."

    set base_cell [molinfo $base_id get {a b c alpha beta gamma}]
    molinfo $merged_id set {a b c alpha beta gamma} $base_cell
    vmdcon -info [format "%-30.30s" "  merged system cell:"] $base_cell

    variable ext [atomselect $merged_id "index >= [$base_sel num]"]
    $ext global
    variable base [atomselect $merged_id "not ([$ext text])"]
    $base global
    set merged [atomselect $merged_id all]
    vmdcon -info "  Base and extension components identified in merged system via '[$base text]' and '[$ext text]'."

    set base_compnames [lsort -unique [$base get $compname]]
    vmdcon -info [format "%-30.30s" "  $compound compunds in base system:"] $base_compnames

    set ext_compnames [lsort -unique [$ext get $compname]]
    vmdcon -info [format "%-30.30s" "  $compound compounds in ext system:"] $ext_compnames

    set merged_compnames [lsort -unique [$merged get $compname]]
    vmdcon -info [format "%-30.30s" "  $compound compounds in merged system:"] $merged_compnames

    vmdcon -info [format "%30s" "#atoms in base system:"] [format "%12d" [$base num]]
    vmdcon -info [format "%30s" "#atoms in ext system:"] [format "%12d" [$ext num]]
    vmdcon -info [format "%30s" "#atoms in merged system:"] [format "%12d" [$merged num]]

    $merged global
    if { ![string equal $autowrap ""]} {
        pbc wrap -molid $merged_id -sel [$merged text] {*}[split $autowrap " "] -all -verbose
    }
    if { ![string equal $autojoin ""]} {
        pbc join {*}[split $autojoin " "] -molid $merged_id -sel [$merged text] -bondlist -all -verbose
    }

    # mol off $base_id
    # mol off $ext_id

    mol rename $merged_id merged

    return $merged_id
}


proc ::MergeTools::overlap { molid base_sel ext_sel } {
    variable overlap_distance
    variable compound

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]

    set overlap [atomselect $molid \
        "([$base_sel text]) and (same $compound as (exwithin $overlap_distance of [$ext_sel text]))"]

    #vmdcon -info [format "%-30.30s" "#atoms in $base_sel overlap with $ext_sel:"] \
    #  [format "%12d" [$overlap num]]

    $overlap global
    return $overlap
}


proc ::MergeTools::report_compounds { molid base_sel {indlvl 0} {linewidth 140} } {
    variable compound
    variable compname

    set base_sel [validate_atomselect $molid $base_sel]

    set immobile_compounds [ mrg -sel $base_sel get immobile ]
    set mobile_compounds [ mrg -sel $base_sel get mobile ]
    set dispensable_compounds [ mrg -sel $base_sel get dispensable ]
    set forbidden_compounds [ mrg -sel $base_sel get forbidden ]


    set complabels {"mobile" "immobile" "dispensable" "forbidden"}
    set complists [ list $mobile_compounds $immobile_compounds $dispensable_compounds $forbidden_compounds ]

    set title [ format "compound report (molid: %d, sel: '%s')" $molid [ $base_sel text ] ]
    set bar [ string repeat "=" [string length $title] ]

    set parstr "- "
    set indstr " "

    set indent [expr 2*($indlvl)]
    set remwidth [expr $linewidth-$indent]
    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

    vmdcon -info ""
    vmdcon -info [format $fmtstr $indstr "$bar"]
    vmdcon -info [format $fmtstr $indstr "$title"]
    vmdcon -info [format $fmtstr $indstr "$bar"]
    vmdcon -info ""

    incr indlvl
    set indent [expr 2*($indlvl)]
    set remwidth [expr $linewidth-$indent]
    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

    vmdcon -info [format $fmtstr $parstr "#atoms ('[ $base_sel text ]'):"] \
        [ format "%12d" [ $base_sel num ] ]

    set unique_compounds [ lsort -unique -integer [ $base_sel get $compound ] ]
    vmdcon -info [format $fmtstr $indstr "#${compound}s ('[ $base_sel text ]'):"] \
        [ format "%12d" [ llength $unique_compounds ] ]

    incr indlvl
    set indent [expr 2*($indlvl)]
    set remwidth [expr $linewidth-$indent]
    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"
    foreach complabel $complabels complist $complists {
        if { [ llength $complist ] == 0 } {
            continue
        }

        if {[$base_sel text] eq "all"} {
            set cur_sel_txt "$compname $complist"
        } else {
            set cur_sel_txt "[$base_sel text] and $compname $complist"
        }

        set cur_sel [atomselect $molid "$cur_sel_txt"]
        #vmdcon -info [format "%2.2s%-138.138s" "- " "$complabel $compname {$complist} in '[$base_sel text]' with '[$ext_sel text]':"]

        vmdcon -info [format $fmtstr $parstr "#atoms in $complabel $compname {$complist} ('[ $cur_sel text ]'):"] \
            [ format "%12d" [ $cur_sel num ] ]

        set unique_compounds [ lsort -unique -integer [ $cur_sel get $compound ] ]
        vmdcon -info [format $fmtstr $indstr "#${compound}s in $complabel $compname {$complist} ('[ $cur_sel text ]'):"] \
            [ format "%12d" [ llength $unique_compounds ] ]

        incr indlvl
        set indent [expr 2*($indlvl)]
        set remwidth [expr $linewidth-$indent]
        set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"
        foreach bc $complist {

            set cur_sel_txt "[$base_sel text] and $compname $bc"
            set cur_sel [atomselect $molid "$cur_sel_txt"]
            vmdcon -info [format $fmtstr $parstr "#atoms in $complabel $compname '$bc' ('[ $cur_sel text ]'):"] \
                [format "%12d" [$cur_sel num]]

            set unique_compounds [ lsort -unique -integer [ $cur_sel get $compound ] ]
            vmdcon -info [format $fmtstr $indstr "#${compound}s in $complabel $compname {$complist} ('[ $cur_sel text ]'):"] \
                [format "%12d" [ llength $unique_compounds ] ]

            unset cur_sel_txt
            unset cur_sel
        }
        incr indlvl -1
    }
    incr indlvl -1
}

proc ::MergeTools::report_overlap { molid base_sel ext_sel } {
    variable overlap_distance
    variable compound
    variable compname

    variable immobile_compounds
    variable mobile_compounds
    variable dispensable_compounds
    variable forbidden_overlap_compounds
    # per default, latter is join of mobile and immobile

    set base_sel [validate_atomselect $molid $base_sel]
    set ext_sel [validate_atomselect $molid $ext_sel]
    # set base_id [validate_atomselect $molid $base_sel]
    # set ext_id [validate_atomselect $molid $ext_sel]
    # TODO: base_id and ext_id must be equal, check

    set complabels {"mobile" "immobile" "dispensable"}
    set complists [list $mobile_compounds $immobile_compounds $dispensable_compounds]

    vmdcon -info ""
    vmdcon -info "==============="
    vmdcon -info "overlap report:"
    vmdcon -info "==============="
    vmdcon -info ""

    # total overlap of base into ext
    set boe [overlap $molid $base_sel $ext_sel]
    vmdcon -info [format "%-140.140s" "#atoms in '[$base_sel text]' overlapping into '[$ext_sel text]':"] \
      [format "%12d" [$boe num]]

    set ext_forbidden_overlap_compound_sel_txt "[$ext_sel text] and $compname $forbidden_overlap_compounds"


    foreach base_complabel $complabels base_complist $complists {
        vmdcon -info [format "%2.2s%-138.138s" "- " "overlap of $base_complabel $compname {$base_complist} in '[$base_sel text]' with '[$ext_sel text]':"]
        foreach bc $base_complist {
            set compound_base_sel_txt "[$base_sel text] and $compname $bc"
            # base compound overlapping ext: bcoe
            set bcoe [overlap $molid $compound_base_sel_txt $ext_sel]
            vmdcon -info [format "%4.4s%-136.136s" "- " "#atoms in '[$base_sel text]' $base_complabel $bc overlapping into '[$ext_sel text]':"] \
                [format "%12d" [$bcoe num]]

            # unique base compound overlapping ext: ubcoe
            set ubcoe [ lsort -unique -integer [ $bcoe get $compound ] ]
            vmdcon -info [format "%4.4s%-136.136s" "  " "#${compound}s in '[$base_sel text]' $base_complabel $bc overlapping into '[$ext_sel text]':"] \
                [format "%12d" [ llength $ubcoe ] ]

            foreach ext_complabel $complabels ext_complist $complists {
                vmdcon -info [format "%6.6s%-134.134s" "- " "overlap of $base_complabel  $compname {$base_complist} in '[$base_sel text]' with $ext_complabel $compname {$ext_complist} in '[$ext_sel text]':"]
                foreach ec $ext_complist {
                    set compound_ext_sel_txt "[$ext_sel text] and $compname $ec"
                    # base compound overlapping ext compound: bcoec
                    set bcoec [overlap $molid $compound_base_sel_txt $compound_ext_sel_txt]
                    vmdcon -info [format "%8.8s%-132.132s" "- " "#atoms in [$base_sel text] $base_complabel $bc overlapping into '[$ext_sel text]' $ext_complabel:"] \
                        [format "%12d" [$bcoec num]]

                    # unique base compound overlapping ext compound: ubcoec
                    set ubcoec [ lsort -unique -integer [ $bcoec get $compound ] ]
                    vmdcon -info [format "%8.8s%-132.132s" "  " "#${compound}s in '[$base_sel text]' $base_complabel $bc overlapping into '[$ext_sel text]' $ext_complabel:"] \
                        [format "%12d" [ llength $ubcoec ] ]

                    unset bcoec
                    unset ubcoec
                }
            }
            unset bcoe
            unset ubcoe
        }
    }
}

proc ::MergeTools::move { molid base_sel ext_sel { position_limits "" } } {
    # position_limits: [ [ float float float ] [ float float float] ]
    #   with the first nexted lists containing the minimum x, y, and z
    #   coordinates and the second containing the corresponding maxima for
    #   randomly generated positions.
    variable maxit

    variable overlap_distance
    variable compound
    variable compname

    variable immobile_compounds
    variable mobile_compounds
    variable dispensable_compounds

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]
    if {$position_limits eq ""} {
        set position_limits [ measure minmax $base_sel ]
    }

    vmdcon -info "Lims: $position_limits"
    vmdcon -info [ format "random position limits: %26s; %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 0 ] ] \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 1 ] ] ]

    set position_origin [ lindex $position_limits 0 ]
    vmdcon -info [ format "random position origin: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$position_origin ] ]
    set position_offset [ vecscale -1.0 [ vecsub {*}$position_limits ] ]
    vmdcon -info [ format "random position offset: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$position_offset ] ]
    # vmd's "residue" and "resid" differ:
    # former is a zero-based consecutive index,
    # latter is identifier as in input file (starts at 1 for standard numbering)
    # but might not be unique

    set forbidden_overlap_compounds [ list {*}$immobile_compounds {*}$mobile_compounds ]

    set nmoved 0
    set ext_forbidden_overlap_compound_sel_txt "[$ext_sel text] and $compname $forbidden_overlap_compounds"
    foreach compnameval $mobile_compounds {
        set base_compound_sel_txt "[$base_sel text] and $compname $compnameval"

        # base compound overlapping ext: bcoe
        set overlapping_compounds [overlap $molid $base_compound_sel_txt $ext_forbidden_overlap_compound_sel_txt]
        vmdcon -info [format "%2.2s%-138.138s" "- " "#atoms in '$base_compound_sel_txt' overlapping with '$ext_forbidden_overlap_compound_sel_txt':"] \
            [format "%12d" [$overlapping_compounds num]]

        # unique base compound overlapping ext: ubcoe
        set unique_overlapping_compound_ids [ lsort -unique -integer [ $overlapping_compounds get $compound ] ]
        vmdcon -info [format "%2.2s%-138.138s" "  " "#${compound}s in '$base_compound_sel_txt' overlapping with '$ext_forbidden_overlap_compound_sel_txt':"] \
            [format "%12d" [ llength $unique_overlapping_compound_ids ] ]

        foreach compid $unique_overlapping_compound_ids {
            vmdcon -info [ format "%4.4s%-136.136s" "- " "teating overlapping $compname $compnameval $compound $compid" ]
            set cur [ atomselect $molid "$compound $compid" ]
            if { [ $cur num ] == 0 } {
                vmdcon -warn [ format "%4.4s%-136.136s" " " "selection empty, already removed!" ]
                continue
            }

            vmdcon -info [ format "%4.4s%-136.136s" " " "#atoms in current $compname $compnameval $compound $compid:" ] \
                [format "%12d" [$cur num]]

            # only try n times to avoid endless loop
            for {set i 0} {$i<$maxit} {incr i} {
                vmdcon -info [format "%6.6s%-134.134s" "- " "iteration $i:" ]

                set cur_overlapping [overlap $molid $cur $ext_forbidden_overlap_compound_sel_txt]
                vmdcon -info [format "%6.6s%-134.134s" " " "#atoms in forbidden overlap with $compname $compnameval $compound $compid:"] \
                    [format "%12d" [$cur_overlapping num]]
                set unique_cur_overlapping_ids [ lsort -unique -integer [ $cur_overlapping get $compound ] ]
                vmdcon -info [format "%6.6s%-134.134s" "  " "#${compound}s in forbidden overlap with $compname $compnameval $compound $compid:"] \
                    [format "%12d" [ llength $unique_cur_overlapping_ids ] ]

                # TODO: outsource mini report
                # report on overlapping molecules
                foreach ext_compnameval $forbidden_overlap_compounds {
                    set ext_compound_sel_txt "[$ext_sel text] and $compname $ext_compnameval"

                    # base compound overlapping ext compound: bcoec
                    set cur_overlapping_compound [overlap $molid $cur $ext_compound_sel_txt]
                    vmdcon -info [format "%8.8s%-132.132s" "- " "#atoms of $compname $ext_compnameval in forbidden overlap with $compname $compnameval $compound $compid:"] \
                        [format "%12d" [$cur_overlapping_compound num]]

                    # unique base compound overlapping ext compound: ubcoec
                    set unique_cur_overlapping_compound_ids [ lsort -unique -integer [ $cur_overlapping_compound get $compound ] ]
                    vmdcon -info [format "%8.8s%-132.132s" "  " "#${compound}s of $compname $ext_compnameval in forbidden overlap with $compname $compnameval $compound $compid:"] \
                        [format "%12d" [ llength $unique_cur_overlapping_compound_ids ] ]

                    unset cur_overlapping_compound
                    unset unique_cur_overlapping_compound_ids
                }

                # check whether position is alright:
                # we allow overlap with dispensable_compounds, which are supposed
                # to be removec subsequently removed, but not with anything else
                # within immobile_compounds and mobile_compounds, i.e. other
                # surfactant chains, counterions or substrate.
                if { [ $cur_overlapping num ] == 0} {
                    vmdcon -info [ format "%4.4s%-136.136s" "- " "found suitable location in iteration $i." ]
                    incr nmoved
                    break
                }

                set random_3vec [ list [ expr rand() ] [ expr rand() ] [ expr rand() ] ]
                vmdcon -info [ format "%6.6s%-134.134s" " " [ format "random 3-vector: %26s" \
                    [ format "%8.4f %8.4f %8.4f" {*}$random_3vec ] ] ]

                set random_position [ vecadd $position_origin [ \
                    vecmul $random_3vec $position_offset] ]
                vmdcon -info [ format "%6.6s%-134.134s" " " [ format "random position: %26s" \
                    [ format "%8.4f %8.4f %8.4f" {*}$random_position ] ] ]

                set cur_position [measure center $cur]
                vmdcon -info [ format "%6.6s%-134.134s" " " [ format "  current position: %26s" \
                    [ format "%8.4f %8.4f %8.4f" {*}$cur_position ] ] ]

                set cur_offset [vecsub $random_position $cur_position]
                vmdcon -info [ format "%6.6s%-134.134s" " " [ format "move by: %26s" \
                    [ format "%8.4f %8.4f %8.4f" {*}$cur_offset ] ] ]
                $cur moveby $cur_offset

                unset cur_overlapping
                unset unique_cur_overlapping_ids
            }
            unset cur
        }
        unset overlapping_compounds
        unset unique_overlapping_compound_ids
    }
    vmdcon -info [ format "successfully moved %3d ${compound}s." $nmoved ]
}

proc ::MergeTools::remove { molid base_sel remove_sel } {
    # creates new molecule from base_sel, with remove_sel removed
    set base_sel [ validate_atomselect $molid $base_sel ]
    set remove_sel [ validate_atomselect $molid $remove_sel ]
    # TODO: check for same molid

    set keep_sel $molid "[$base_sel text] and not [$remove_sel text]"
    variable keep_id [::TopoTools::selections2mol $keep_sel]

    # copy periodic box
    molinfo $keep_id set {a b c alpha beta gamma} \
        [molinfo $molid get {a b c alpha beta gamma}]
}

interp alias {} mrg {} ::MergeTools::mrg
package provide mergetools $::MergeTools::version
