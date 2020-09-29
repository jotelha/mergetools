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
# Refer to README.md for more iformation and examples.
namespace eval ::MergeTools:: {
    variable version 0.1

    package require topotools
    package require pbctools
    package require struct::set

    # distance within which molecules are regarded as overlapping
    variable overlap_distance 2.0
    variable maxit 1000

    variable autowrap "-nocompound"
    variable autojoin "residue"
    variable bondlist "-bondlist"

    variable compound "residue"
    variable compname "resname"

    variable negative_keywords {"no" "none" "off" "false" "0"}
    variable positive_keywords {"yes" "on" "true" "1"}


    variable valid_autowrap_vals [list \
        {*}$negative_keywords \
        "atom" "res" "resid" "residue" "seg" "segid" "chain" "fragment"]

    variable valid_autojoin_vals [list \
        {*}$negative_keywords \
        "seg" "segid" "res" "resid" "residue" \
        "chain" "fragment" "bonded" "connected"]

    variable valid_compound_vals [list \
        "res" "resid" "residue" "seg" "segid" "chain" "fragment"]

    variable valid_compname_vals [list \
        "resname" "segname" ]
    # what else?

    variable immobile_compounds {}
    variable mobile_compounds -1
    variable dispensable_compounds {SOL}
    variable forbidden_compounds -1

    variable base
    variable ext

    variable skip_zero 1
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
    vmdcon -info "        Please refer to the README.md file shipped with this"
    vmdcon -info "        package for detailed usage information."
    vmdcon -info ""
    vmdcon -info "common flags (not implemented):"
    vmdcon -info ""
    vmdcon -info "  -molid     <num>|top    molecule id (default: 'top')"
    vmdcon -info "  -sel       <selection>  atom selection function or text (default: 'all')"
    vmdcon -info ""
    vmdcon -info "sample usage:"
    vmdcon -info ""
    vmdcon -info "  vmd> package require mergetool"
    vmdcon -info "  vmd>"
    vmdcon -info "  vmd> set substrate_system_id \[mol new substrate.gro\]"
    vmdcon -info "  vmd> set indenter_system_id \[mol new indenter.gro\]"
    vmdcon -info "  vmd> set substrate_system \[atomselect \$substrate_system_id all\]"
    vmdcon -info "  vmd> set indenter_system \[atomselect \$indenter_system_id all\]"
    vmdcon -info "  vmd>"
    vmdcon -info "  vmd> mrg set autowrap residue autojoin off"
    vmdcon -info "  vmd> mrg set dispensable SOL"
    vmdcon -info "  vmd> mrg set immobile AUM"
    vmdcon -info "  vmd>"
    vmdcon -info "  vmd> set merged_id \[mrg -sel \$substrate_system merge \$indenter_system\]"
    vmdcon -info "  vmd> set base \[mrg get base\]"
    vmdcon -info "  vmd> set ext \[mrg get ext\]"
    vmdcon -info "  vmd> set merged \[atomselect \$merged_id all\]"
    vmdcon -info "  vmd>"
    vmdcon -info "  vmd> set overlap_to_move \[mrg -sel \$base overlap mobile \$ext]"
    vmdcon -info "  vmd> mrg -sel \$merged move \$overlap_to_move"
    vmdcon -info "  vmd>"
    vmdcon -info "  vmd> set overlap_to_remove \[mrg -sel \$base overlap removable \$ext]"
    vmdcon -info "  vmd> set nonoverlap_id \[mrg -sel \$merged remove \$overlap_to_remove\]"
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
    variable positive_keywords
    variable valid_autowrap_vals
    variable valid_autojoin_vals
    variable valid_compound_vals
    variable valid_compname_vals

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
      "remove"
      "subtract"
      "get"
      "set"
      "report"
      "help"
    }

      # valid properties.
      set validprop {
        "autojoin"
        "autowrap"
        "base"
        "compound"
        "compname"
        "compnames"
        "dispensable"
        "ext"
        "forbidden"
        "immobile"
        "mobile"
        "maxit"
        "overlap"
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

                bondlist {
                    variable bondlist
                    if {[string equal $bondlist ""]} {
                        set retval 0
                    } else {
                        set retval 1
                    }
                }

                compound {
                    variable compound
                    set retval $compound
                }

                compname {
                    variable compname
                    set retval $compname
                }

                overlap {
                    variable overlap_distance
                    set retval $overlap_distance
                }

                maxit {
                    variable maxit
                    set retval $maxit
                }

                compnames {
                    set retval [ compnames $molid $sel ]
                }

                base {
                    variable base
                    set retval $base
                    $retval uplevel 1
                }

                ext {
                    variable ext
                    set retval $ext
                    $retval uplevel 1
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

                    bondlist {
                        set val [lindex $newargs 0]
                        if {$val in $positive_keywords} {
                            variable bondlist "-bondlist"
                        } else {
                            variable bondlist ""
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

                    compname {
                        set val [lindex $newargs 0]
                        if {$val ni $valid_compname_vals} {
                        		vmdcon -err "'merge set $key' allows choice from '$valid_compname_vals'"
                            return
                        }
                        variable compname $val
                        set newargs [lrange $newargs 1 end]
                    }

                    overlap {
                        set val [lindex $newargs 0]
                        variable overlap_distance $val
                        set newargs [lrange $newargs 1 end]
                    }

                    maxit {
                      set val [lindex $newargs 0]
                      variable maxit $val
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
                            set valid_vals [ compnames $molid $sel ]
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
            set key [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            switch -nocase -- $key {
                removable {
                    set retval [overlap_removable $molid $sel [lindex $newargs 0]]
                }
                dispensable {
                    set retval [overlap_dispensable $molid $sel [lindex $newargs 0]]
                }
                forbidden {
                    set retval [overlap_forbidden $molid $sel [lindex $newargs 0]]
                }
                mobile {
                    set retval [overlap_mobile $molid $sel [lindex $newargs 0]]
                }
                twoway -
                bi -
                bidirectional {
                    set retval [overlap_bidirectional $molid $sel [lindex $newargs 0]]
                }
                oneway -
                uni -
                unidirectional {
                    set retval [overlap $molid $sel [lindex $newargs 0]]
                }
                internal {
                    set retval [overlap_internal $molid $sel]
                }
                default {
                    set retval [overlap $molid $sel $key]
                }
            }
            # we have to manualy elevate atomelect procs
            $retval uplevel 1
        }
        "move" {
            set ext [lindex $newargs 0]
            set newargs [lrange $newargs 1 end]
            if {[llength $newargs] == 0} {
                set retval [move $molid $sel $ext]
            } else {
                # position_limits explicitly specified
                set position_limits [lindex $newargs 0]
                set retval [move $molid $sel $ext $position_limits]
            }
        }
        "subtract" {
            set retval [subtract $molid $sel [lindex $newargs 0]]
        }
        "remove" {
            set retval [remove $molid $sel [lindex $newargs 0]]
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
                    report_overlap $molid $sel [lindex $newargs 0]
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

    # for selecting overlaps, it is desirable to have all atoms within the cell
    if { ![string equal $autowrap ""]} {
        pbc wrap -molid $merged_id -sel [$merged text] {*}[split $autowrap " "] -all -verbose
    }

    $merged delete
    return $merged_id
}


proc ::MergeTools::overlap { molid base_sel ext_sel {ex "ex"} } {
    # select compounds in ext_sel that overlap with compounds in base_sel
    variable overlap_distance
    variable compound

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]

    set overlap [ atomselect $molid \
        "([$ext_sel text]) and (same $compound as (${ex}within $overlap_distance of [$base_sel text]))" ]

    # TODO: nasty legacy behavior in VMD, see https://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/21856.html
    # $overlap global won't work for the caller, thus upproc 1 necessary.
    $overlap uplevel 1
    return $overlap
}


proc ::MergeTools::overlap_removable { molid base_sel ext_sel {ex "ex"} } {
    # finds dispensable overlap whose removal will render all remaining
    # dispensable compounds in the system non-overlapping
    # strategy as follows: first find all dispensable compounds in ext that
    # overlap with any compound in base and mark them for removal; next, find
    # all dispensable compounds in base that overlap with any compound in ext
    # not marked for removal; return union of both
    variable compound
    variable compname

    set dispensable_compounds [ mrg -sel $base_sel get dispensable ]

    set ext_dispensable [atomselect $molid "$compname $dispensable_compounds and ([$ext_sel text])"]
    set base_dispensable [atomselect $molid "$compname $dispensable_compounds and ([$base_sel text])"]

    set ext_to_remove [overlap $molid $base_sel $ext_dispensable]
    set ext_to_keep [atomselect $molid "([$ext_sel text]) and not ([$ext_to_remove text])"]

    set base_to_remove [overlap $molid $ext_to_keep $base_dispensable]

    set overlap [atomselect $molid "([$base_to_remove text]) or ([$ext_to_remove text])"]

    $ext_dispensable delete
    $base_dispensable delete
    $ext_to_remove delete
    $ext_to_keep delete
    $base_to_remove delete

    $overlap uplevel 1
    return $overlap
}


proc ::MergeTools::overlap_mobile { molid base_sel ext_sel {ex "ex"} } {
    # select union of mobile compounds in ext_sel that overlap with forbidden
    # compounds in base_sel and mobile compounds in base_sel that overlap with
    # forbidden compounds in ext_sel
    variable overlap_distance
    variable compname

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]

    set mobile_compounds [ mrg -sel $base_sel get mobile ]
    set forbidden_compounds [ mrg -sel $base_sel get forbidden ]

    set base_forbidden [ atomselect $molid "([$base_sel text]) and ($compname $forbidden_compounds)" ]
    set ext_forbidden [ atomselect $molid "([$ext_sel text]) and ($compname $forbidden_compounds)" ]
    set base_mobile [ atomselect $molid "([$base_sel text]) and ($compname $mobile_compounds)" ]
    set ext_mobile [ atomselect $molid "([$ext_sel text]) and ($compname $mobile_compounds)" ]

    # backward overlap: from ext into base
    set backward_overlap [ overlap $molid $base_forbidden $ext_mobile ]
    # forward overlap: from base into ext
    set forward_overlap [ overlap $molid $ext_forbidden $base_mobile ]

    set overlap [ atomselect $molid "([ $backward_overlap text ]) or ([ $forward_overlap text ])" ]

    $base_forbidden delete
    $ext_forbidden delete
    $base_mobile delete
    $ext_mobile delete

    $backward_overlap delete
    $forward_overlap delete

    $overlap uplevel 1
    return $overlap
}


proc ::MergeTools::overlap_forbidden { molid base_sel ext_sel {ex "ex"} } {
    # select union of any compound in ext_sel that overlap with forbidden
    # compounds in base_sel and any compound in base_sel that overlap with
    # forbidden compounds in ext_sel
    variable overlap_distance
    variable compname

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]

    set forbidden_compounds [ mrg -sel $base_sel get forbidden ]

    set base_forbidden [ atomselect $molid "([$base_sel text]) and ($compname $forbidden_compounds)" ]
    set ext_forbidden [ atomselect $molid "([$ext_sel text]) and ($compname $forbidden_compounds)" ]

    # backward overlap: from ext into base
    set backward_overlap [ overlap $molid $base_forbidden $ext_sel ]
    # forward overlap: from base into ext
    set forward_overlap [ overlap $molid $ext_forbidden $base_sel ]

    set overlap [ atomselect $molid "([ $backward_overlap text ]) or ([ $forward_overlap text ])" ]

    $base_forbidden delete
    $ext_forbidden delete
    $forward_overlap delete
    $backward_overlap delete
    $overlap uplevel 1
    return $overlap
}


proc ::MergeTools::overlap_dispensable { molid base_sel ext_sel {ex "ex"} } {
    # select union of dispensable compound in ext_sel that overlap with any
    # compounds in base_sel and dispensable compound in base_sel that overlap
    # with any compounds in ext_sel
    variable overlap_distance
    variable compname

    set base_sel [ validate_atomselect $molid $base_sel ]
    set ext_sel [ validate_atomselect $molid $ext_sel ]

    set dispensable_compounds [ mrg -sel $base_sel get dispensable ]

    set base_dispensable [ atomselect $molid "([$base_sel text]) and ($compname $dispensable_compounds)" ]
    set ext_dispensable [ atomselect $molid "([$ext_sel text]) and ($compname $dispensable_compounds)" ]

    # backward overlap: from ext into base
    set backward_overlap [ overlap $molid $base_sel $ext_dispensable ]
    # forward overlap: from base into ext
    set forward_overlap [ overlap $molid $ext_sel $base_dispensable ]

    set overlap [ atomselect $molid "([ $backward_overlap text ]) or ([ $forward_overlap text ])" ]

    $forward_overlap delete
    $backward_overlap delete
    $overlap uplevel 1
    return $overlap
}

proc ::MergeTools::overlap_bidirectional { molid base_sel ext_sel {ex "ex"} } {
    # select union of compounds in ext_sel that overlap with compounds in
    # base_sel and compounds in base_sel that overlap with compounds in ext_sel
    variable overlap_distance
    variable compound

    # backward overlap: from ext into base
    set backward_overlap [ overlap $molid $base_sel $ext_sel ]
    # forward overlap: from base into ext
    set forward_overlap [ overlap $molid $ext_sel $base_sel ]

    set overlap [ atomselect $molid "([ $backward_overlap text ]) or ([ $forward_overlap text ])" ]

    $forward_overlap delete
    $backward_overlap delete
    $overlap uplevel 1
    return $overlap
}


proc ::MergeTools::overlap_internal { molid base_sel {indlvl 0} {linewidth 180} } {
    # finds a subselection whose removal will render base_sel overlap-free
    variable compound
    variable compname

    set base_sel [ validate_atomselect $molid $base_sel ]

    set parstr "- "
    set indstr " "

    set indent [expr 2*($indlvl)]
    set remwidth [expr $linewidth-$indent]
    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

    incr indlvl
    set max_penalty 0
    set worst_compid -1

    set unique_overlap_compound_ids {}
    set niter 0
    incr indlvl
    while 1 {
        set worst_compid -1
        # find compound with worst overlap
        if {[$base_sel text] eq "all"} {
            if { [llength $unique_overlap_compound_ids] > 0 } {
                set overlap_sel_txt "$compname $unique_overlap_compound_ids"
                set nonoverlap_sel_txt "not $compname $unique_overlap_compound_ids"
            } else {
                set overlap_sel_txt "none"
                set nonoverlap_sel_txt "all"
            }
        } else {
            if { [llength $unique_overlap_compound_ids] > 0 } {
                set overlap_sel_txt "$compname $unique_overlap_compound_ids and ([$base_sel text])"
                set nonoverlap_sel_txt "not $compname $unique_overlap_compound_ids and ([$base_sel text])"
            } else {
                set overlap_sel_txt "none"
                set nonoverlap_sel_txt "([$base_sel text])"
            }
        }

        set overlap_sel [ atomselect $molid "$overlap_sel_txt" ]
        set nonoverlap_sel [ atomselect $molid "$nonoverlap_sel_txt" ]

        set unique_compound_ids [ lsort -unique -integer [ $nonoverlap_sel get $compound ] ]
        set tot_count [llength $unique_compound_ids]
        set count 0
        incr indlvl
        vmdcon -nonewline -info [format "%${indent}.${indent}s#iter $niter progress:" $indstr]
        flush stdout
        foreach compid $unique_compound_ids {
            if {[$nonoverlap_sel text] eq "all"} {
                set cur_sel_txt "$compname $compid"
                set cur_neg_sel_txt "not $compname $compid"
            } else {
                set cur_sel_txt "$compname $compid and ([$base_sel text])"
                set cur_neg_sel_txt "not $compname $compid and ([$nonoverlap_sel text])"
            }

            set cur_sel [ atomselect $molid "$cur_sel_txt" ]
            set cur_neg_sel [ atomselect $molid "$cur_neg_sel_txt" ]

            set cur_overlap [ overlap $molid $cur_sel $cur_neg_sel ]

            set unique_cur_overlap_compound_ids [ lsort -unique -integer [ $cur_overlap get $compound ] ]
            set cur_penalty [llength $unique_cur_overlap_compound_ids]
            if { $cur_penalty > $max_penalty } {
                set worst_compid $compid
                set worst_sel $cur_sel
                set max_penalty $cur_penalty
            }

            $cur_sel delete
            $cur_neg_sel delete
            $cur_overlap delete

            puts -nonewline "."
            flush stdout
            incr count
        }
        incr indlvl -1
        puts ""

        if { $worst_compid == -1 } {
            vmdcon -info [ format $fmtstr $indstr "Iteration $iter: remaining selection overlap-free." ]
            $overlap_sel delete
        } else {
            set worst_compnameval [ lsort -unique [$worst_sel get $compname] ]
            vmdcon -info [ format $fmtstr $indtsr "Iteration $niter: identified worst overlap of $worst_penalty at $compname $worst_compnameval $compound $worst_compid." ]
            lappend unique_overlap_compound_ids $worst_compid
            $overlap_sel delete
            $nonoverlap_sel delete
        }
        incr niter
    }
    incr indlvl -1

    $overlap_sel uplevel 1
    return $overlap_sel
}


proc ::MergeTools::report_compounds { molid base_sel {indlvl 0} {linewidth 180} } {
    variable compound
    variable compname
    variable skip_zero

    set loc_base_sel [validate_atomselect $molid $base_sel]

    set immobile_compounds [ mrg -sel $loc_base_sel get immobile ]
    set mobile_compounds [ mrg -sel $loc_base_sel get mobile ]
    set dispensable_compounds [ mrg -sel $loc_base_sel get dispensable ]
    set forbidden_compounds [ mrg -sel $loc_base_sel get forbidden ]

    set complabels {"mobile" "immobile" "dispensable" "forbidden"}
    set complists [ list $mobile_compounds $immobile_compounds $dispensable_compounds $forbidden_compounds ]

    set title [ format "compound report (molid: %d, sel: '%s')" $molid [ $loc_base_sel text ] ]
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

    vmdcon -info [format $fmtstr $parstr "#atoms ('[ $loc_base_sel text ]'):"] \
        [ format "%12d" [ $loc_base_sel num ] ]

    set unique_compounds [ lsort -unique -integer [ $loc_base_sel get $compound ] ]
    vmdcon -info [format $fmtstr $indstr "#${compound}s ('[ $loc_base_sel text ]'):"] \
        [ format "%12d" [ llength $unique_compounds ] ]

    incr indlvl
    foreach complabel $complabels complist $complists {
        set indent [expr 2*($indlvl)]
        set remwidth [expr $linewidth-$indent]
        set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

        if { [ llength $complist ] == 0 } {
            continue
        }

        if {[$loc_base_sel text] eq "all"} {
            set cur_sel_txt "$compname $complist"
        } else {
            set cur_sel_txt "[$loc_base_sel text] and $compname $complist"
        }

        set cur_sel [atomselect $molid "$cur_sel_txt"]
        if { $skip_zero && [$cur_sel num] == 0 } {
            continue
        }
        set unique_compounds [ lsort -unique -integer [ $cur_sel get $compound ] ]

        vmdcon -info [format $fmtstr $parstr "#atoms in $complabel $compname {$complist} ('[ $cur_sel text ]'):"] \
            [ format "%12d" [ $cur_sel num ] ]
        vmdcon -info [format $fmtstr $indstr "#${compound}s in $complabel $compname {$complist} ('[ $cur_sel text ]'):"] \
            [ format "%12d" [ llength $unique_compounds ] ]
        $cur_sel delete

        incr indlvl
        foreach bc $complist {
            set indent [expr 2*($indlvl)]
            set remwidth [expr $linewidth-$indent]
            set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

            set cur_sel_txt "[$loc_base_sel text] and $compname $bc"
            set cur_sel [atomselect $molid "$cur_sel_txt"]
            if { $skip_zero && [$cur_sel num] == 0 } {
                continue
            }
            set unique_compounds [ lsort -unique -integer [ $cur_sel get $compound ] ]

            vmdcon -info [format $fmtstr $parstr "#atoms in $complabel $compname '$bc' ('[ $cur_sel text ]'):"] \
                [format "%12d" [$cur_sel num]]
            vmdcon -info [format $fmtstr $indstr "#${compound}s in $complabel $compname '$bc' ('[ $cur_sel text ]'):"] \
                [format "%12d" [ llength $unique_compounds ] ]


            $cur_sel delete
        }
        incr indlvl -1
    }
    incr indlvl -1
}


proc ::MergeTools::report_overlap { molid base_sel ext_sel {indlvl 0} {linewidth 180} } {
    # report compounds in ext_sel that overlap with compounds in base_sel
    variable overlap_distance
    variable compound
    variable compname
    variable skip_zero

    set immobile_compounds [mrg get immobile]
    set mobile_compounds [mrg get mobile]
    set dispensable_compounds [mrg get dispensable]
    set forbidden_compounds [mrg get forbidden]
    # per default, latter is join of mobile and immobile

    set base_sel [validate_atomselect $molid $base_sel]
    set ext_sel [validate_atomselect $molid $ext_sel]
    # TODO: base_id and ext_id must be equal, check

    set complabels {"mobile" "immobile" "dispensable" "forbidden"}
    set complists [list $mobile_compounds $immobile_compounds $dispensable_compounds $forbidden_compounds]

    set title [ format "overlap report (molid: %d, '%s' in '%s')" $molid [ $ext_sel text ] [ $base_sel text ] ]
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

    # eob: ext overlap with  base
    set eob [overlap $molid $base_sel $ext_sel]
    # ueob: unique compounds in ext overlap with base
    set ueob [ lsort -unique -integer [ $eob get $compound ] ]
    vmdcon -info [ format $fmtstr $parstr "#atoms overlapping:" ] [ format "%12d" [ $eob num ] ]
    vmdcon -info [ format $fmtstr $indstr "#${compound}s overlapping:" ] [ format "%12d" [ llength $ueob ] ]

    incr indlvl
    foreach ext_complabel $complabels ext_complist $complists {
        set indent [expr 2*($indlvl)]
        set remwidth [expr $linewidth-$indent]
        set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

        if { [ llength $ext_complist ] == 0 } {
            continue
        }
        if {[$ext_sel text] eq "all"} {
            set ext_complist_sel_txt "$compname $ext_complist"
        } else {
            set ext_complist_sel_txt "[$ext_sel text] and $compname $ext_complist"
        }
        # elob: ext compname list overlap with base
        set elob [overlap $molid $base_sel $ext_complist_sel_txt]
        if { $skip_zero && [$elob num] == 0 } {
            continue
        }
        set uelob [ lsort -unique -integer [ $elob get $compound ] ]

        vmdcon -info [ format $fmtstr $parstr "#atoms in $ext_complabel $compname {$ext_complist} overlapping:" ] [ format "%12d" [ $elob num ] ]
        vmdcon -info [ format $fmtstr $indstr "#${compound}s in $ext_complabel $compname {$ext_complist} overlapping:" ] [ format "%12d" [ llength $uelob ] ]

        incr indlvl
        foreach ec $ext_complist {
            set indent [expr 2*($indlvl)]
            set remwidth [expr $linewidth-$indent]
            set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

            if {[$ext_sel text] eq "all"} {
                set ext_compname_sel_txt "$compname $ec"
            } else {
                set ext_compname_sel_txt "[$ext_sel text] and $compname $ec"
            }
            # ecob: ext compname overlap with base
            set ecob [ overlap $molid $base_sel $ext_compname_sel_txt ]
            if { $skip_zero && [$ecob num] == 0 } {
                continue
            }
            # uecob: unique compounds inext compname overlap with base
            set uecob [ lsort -unique -integer [ $ecob get $compound ] ]
            vmdcon -info [ format $fmtstr $parstr "#atoms in $ext_complabel $compname '$ec' overlapping:" ] [ format "%12d" [ $ecob num ] ]
            vmdcon -info [ format $fmtstr $indstr "#${compound}s in $ext_complabel $compname '$ec' overlapping:" ] [ format "%12d" [ llength $uecob ] ]

            incr indlvl
            foreach base_complabel $complabels base_complist $complists {
                set indent [expr 2*($indlvl)]
                set remwidth [expr $linewidth-$indent]
                set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

                if { [ llength $base_complist ] == 0 } {
                    continue
                }
                if {[$base_sel text] eq "all"} {
                    set base_complist_sel_txt "$compname $base_complist"
                } else {
                    set base_complist_sel_txt "[$base_sel text] and $compname $base_complist"
                }
                # ecobl: ext compname overlap with base compname list
                set ecobl [ overlap $molid $base_complist_sel_txt $ext_compname_sel_txt ]
                if { $skip_zero && [$ecobl num] == 0 } {
                    continue
                }
                # elob: unique compname in ext overlap with base compname list
                set uecobl [ lsort -unique -integer [ $ecobl get $compound ] ]

                vmdcon -info [ format $fmtstr $parstr "#atoms in $ext_complabel $compname '$ec' overlapping $base_complabel $compname {$base_complist}:" ] \
                    [ format "%12d" [ $ecobl num ] ]
                vmdcon -info [ format $fmtstr $indstr "#${compound}s in $ext_complabel $compname '$ec' overlapping $base_complabel $compname {$base_complist}:" ] \
                    [ format "%12d" [ llength $uecobl ] ]

                incr indlvl
                foreach bc $base_complist {
                    set indent [expr 2*($indlvl)]
                    set remwidth [expr $linewidth-$indent]
                    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

                    if {[$base_sel text] eq "all"} {
                        set base_compname_sel_txt "$compname $bc"
                    } else {
                        set base_compname_sel_txt "[$base_sel text] and $compname $bc"
                    }
                    # ecobc: ext compound overlap with base compound
                    set ecobc [overlap $molid $base_compname_sel_txt $ext_compname_sel_txt]
                    if { $skip_zero && [$ecobc num] == 0 } {
                        continue
                    }
                    # uecobc: unique compname in  ext compound overlapping with base compound
                    set uecobc [ lsort -unique -integer [ $ecobc get $compound ] ]

                    vmdcon -info [ format $fmtstr $parstr "#atoms in $ext_complabel $compname '$ec' overlapping $base_complabel $compname '$bc':" ] \
                        [ format "%12d" [ $ecobc num ] ]
                    vmdcon -info [ format $fmtstr $indstr "#${compound}s in $ext_complabel $compname '$ec' overlapping $base_complabel $compname '$bc':" ] \
                        [ format "%12d" [ llength $uecobc ] ]

                    $ecobc delete
                    $uecobc delete
                }
                incr indlvl -1
                $ecobl delete
                $uecobl delete
            }
            incr indlvl -1
            $ecob delete
            $uecob delete
        }
        incr indlvl -1
        $elob delete
        $uelob delete
    }
    incr indlvl -1
    $eob delete
    $ueob delete
}


proc ::MergeTools::move { molid base_sel overlap_sel { position_limits "" } {indlvl 0} {linewidth 180} } {
    # position_limits: [ [ float float float ] [ float float float] ]
    #   with the first nexted lists containing the minimum x, y, and z
    #   coordinates and the second containing the corresponding maxima for
    #   randomly generated positions.
    # variable maxit

    variable overlap_distance
    variable compound
    variable compname
    variable autojoin
    variable autowrap
    variable bondlist

    set base_sel [ validate_atomselect $molid $base_sel ]
    set overlap_sel [ validate_atomselect $molid $overlap_sel ]

    # when moving around molecules, it is desirable to have them joint
    if { ![string equal $autojoin ""]} {
        pbc join {*}[split $autojoin " "] -molid $molid -sel [$overlap_sel text] $bondlist -all -verbose
    }

    set mobile_compounds [ mrg -sel $base_sel get mobile ]
    set forbidden_compounds [ mrg -sel $base_sel get forbidden ]

    if {$position_limits eq ""} {
        set position_limits [ measure minmax $base_sel ]
    }

    set parstr "- "
    set indstr " "

    set indent [expr 2*($indlvl)]
    set remwidth [expr $linewidth-$indent]
    set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

    vmdcon -info [ format $fmtstr $indstr "Lims: $position_limits" ]
    vmdcon -info [ format $fmtstr $indstr [format "random position limits: %26s; %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 0 ] ] \
        [ format "%8.4f %8.4f %8.4f" {*}[ lindex $position_limits 1 ] ] ] ]

    set position_origin [ lindex $position_limits 0 ]
    vmdcon -info [ format $fmtstr $indstr [ format "random position origin: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$position_origin ] ] ]
    set position_offset [ vecscale -1.0 [ vecsub {*}$position_limits ] ]
    vmdcon -info [ format $fmtstr $indstr [ format "random position offset: %26s" \
        [ format "%8.4f %8.4f %8.4f" {*}$position_offset ] ] ]
    # vmd's "residue" and "resid" differ:
    # former is a zero-based consecutive index,
    # latter is identifier as in input file (starts at 1 for standard numbering)
    # but might not be unique

    set nmoved 0
    set nfailed 0

    if {[$base_sel text] eq "all"} {
        set base_forbidden_sel [atomselect $molid "$compname $forbidden_compounds"]
    } else {
        set base_forbidden_sel [atomselect $molid "[$base_sel text] and $compname $forbidden_compounds"]
    }
    vmdcon -info ""
    vmdcon -info [ format $fmtstr $indstr "forbidden ${compound}s" ]
    vmdcon -info ""
    report_compounds $molid $base_forbidden_sel $indlvl

    # explicit loop here meant to ensure specific order in movements
    incr indlvl
    foreach compnameval $mobile_compounds {
        set indent [expr 2*($indlvl)]
        set remwidth [expr $linewidth-$indent]
        set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

        set overlap_compname_sel [atomselect $molid "[$overlap_sel text] and $compname $compnameval" ]

        set unique_overlap_compound_ids [ lsort -unique -integer [ $overlap_compname_sel get $compound ] ]

        incr indlvl
        foreach compid $unique_overlap_compound_ids {
            set indent [expr 2*($indlvl)]
            set remwidth [expr $linewidth-$indent]
            set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

            vmdcon -info [ format $fmtstr $parstr "$compname $compnameval $compound $compid" ]

            set cur_move_sel [ atomselect $molid "$compound $compid" ]
            set cur_forbidden_sel [ atomselect $molid "[$base_forbidden_sel text] and not $compound $compid" ]

            if { [ $cur_move_sel num ] == 0 } {
                vmdcon -error [ format $fmtstr $indstr "selection empty, already removed!" ]
                break
            }

            set ret [ move_one $molid $cur_forbidden_sel $cur_move_sel $position_origin $position_offset $indlvl $linewidth]
            if { $ret } { incr nmoved } else { incr nfailed }

            $cur_move_sel delete
            $cur_forbidden_sel delete
        }
        incr indlvl -1

        $overlap_compname_sel delete
    }
    incr indlvl -1
    if { $nfailed > 0 } {
        vmdcon -warn [ format $fmtstr $indstr [ format "faile moving %3d ${compound}s." $nfailed ] ]
    }
    vmdcon -info [ format $fmtstr $indstr [ format "successfully moved %3d ${compound}s." $nmoved ] ]

    # TODO: it might be necessary to wrap after each iteration
    if { ![string equal $autowrap ""]} {
        pbc wrap -molid $molid -sel [$overlap_sel text] {*}[split $autowrap " "] -all -verbose
    }
}


proc ::MergeTools::move_one { molid base_sel move_sel position_origin position_offset {indlvl 0} {linewidth 180} } {
    # move move_sel until it has no overlap with base_sel anymore
    variable maxit

    set parstr "- "
    set indstr " "

    incr indlvl
    # only try n times to avoid endless loop
    for {set i 0} {$i<$maxit} {incr i} {
        set indent [expr 2*($indlvl)]
        set remwidth [expr $linewidth-$indent]
        set fmtstr "%${indent}.${indent}s%-${remwidth}.${remwidth}s"

        vmdcon -info [ format $fmtstr $parstr "iteration $i:" ]
        # report_overlap $molid $base_sel $move_sel $indlvl $linewidth

        set cur_overlapping [ overlap $molid $base_sel $move_sel ]

        # check whether position is alright:
        # we allow overlap with dispensable_compounds, which are supposed
        # to be subsequently removed, but not with anything else
        # within immobile_compounds and mobile_compounds, i.e. other
        # surfactant chains, counterions or substrate.
        if { [ $cur_overlapping num ] == 0} {
            vmdcon -info [ format $fmtstr $indstr "found suitable location  in iteration $i." ]
            $cur_overlapping delete
            return 1
        }

        set random_3vec [ list [ expr rand() ] [ expr rand() ] [ expr rand() ] ]
        vmdcon -info [ format $fmtstr $indstr [ format "random 3-vector: %26s" \
            [ format "%8.4f %8.4f %8.4f" {*}$random_3vec ] ] ]

        set random_position [ vecadd $position_origin [ \
            vecmul $random_3vec $position_offset] ]
        vmdcon -info [ format $fmtstr $indstr [ format "random position: %26s" \
            [ format "%8.4f %8.4f %8.4f" {*}$random_position ] ] ]

        set cur_position [measure center $move_sel]
        vmdcon -info [ format $fmtstr $indstr [ format "current position: %26s" \
            [ format "%8.4f %8.4f %8.4f" {*}$cur_position ] ] ]

        set cur_offset [vecsub $random_position $cur_position]
        vmdcon -info [ format $fmtstr $indstr [ format "move by: %26s" \
            [ format "%8.4f %8.4f %8.4f" {*}$cur_offset ] ] ]

        $move_sel moveby $cur_offset

        $cur_overlapping delete
    }
    incr indlvl -1
    return 0
}

proc ::MergeTools::subtract { molid base_sel remove_sel } {
    # creates new molecule from base_sel, with remove_sel removed
    variable compname
    variable compound
    variable autojoin
    set base_sel [ validate_atomselect $molid $base_sel ]
    set remove_sel [ validate_atomselect $molid $remove_sel ]
    # TODO: check for same molid

    # weird behavior in directly using negated selection text here, i.e.
    # set keep_sel [atomselect $molid "[$base_sel text] and not [$remove_sel text]"]
    # my susppect is the exwithin mechanism
    # TODO: check source of issue and adapt selection texts elsewhere in this
    # package if necessary. Workaround for now: Selections based on compound.
    set base_compounds [lsort -unique -integer [$base_sel get $compound]]
    set remove_compounds [lsort -unique -integer [$remove_sel get $compound]]

    set keep_compounds [struct::set difference $base_compounds $remove_compounds]
    set keep_sel [atomselect $molid "$compound $keep_compounds"]

    $keep_sel uplevel 1
    return $keep_sel
}


proc ::MergeTools::remove { molid base_sel remove_sel } {
    # creates new molecule from base_sel, with remove_sel removed
    variable compname
    variable compound
    variable autojoin
    set base_sel [ validate_atomselect $molid $base_sel ]
    set remove_sel [ validate_atomselect $molid $remove_sel ]
    # TODO: check for same molid

    # weird behavior in directly using negated selection text here, i.e.
    # set keep_sel [atomselect $molid "[$base_sel text] and not [$remove_sel text]"]
    # my susppect is the exwithin mechanism
    # TODO: check source of issue and adapt selection texts elsewhere in this
    # package if necessary. Workaround for now: Selections based on compound.
    set base_compounds [lsort -unique -integer [$base_sel get $compound]]
    set remove_compounds [lsort -unique -integer [$remove_sel get $compound]]

    set keep_compounds [lmap item $base_compounds {
        if {$item ni $remove_compounds} {set item} else {continue}
    }]
    # set keep_compounds [struct::set difference $base_compounds $remove_compounds]
    set keep_sel [atomselect $molid "$compound $keep_compounds"]
    set keep_id [::TopoTools::selections2mol [list $keep_sel]]

    # copy periodic box
    molinfo $keep_id set {a b c alpha beta gamma} \
        [molinfo $molid get {a b c alpha beta gamma}]

    # it might be desired to join split residues as this is usually the last
    # step of mergin systems
    # TODO: autowrapping and joining might be desirable at other places as well
    if { ![string equal $autojoin ""]} {
        pbc join {*}[split $autojoin " "] -molid $keep_id -sel all -bondlist -all -verbose
    }

    $keep_sel delete
    return $keep_id
}


proc ::MergeTools::lmap args {
    set body [lindex $args end]
    set args [lrange $args 0 end-1]
    set n 0
    set pairs [list]
    foreach {varnames listval} $args {
        set varlist [list]
        foreach varname $varnames {
            upvar 1 $varname var$n
            lappend varlist var$n
            incr n
        }
        lappend pairs $varlist $listval
    }
    set temp [list]
    foreach {*}$pairs {
        lappend temp [uplevel 1 $body]
    }
    set temp
}


interp alias {} mrg {} ::MergeTools::mrg
package provide mergetools $::MergeTools::version
