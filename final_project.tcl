#############################################################################
# Generated by PAGE version 8.0
#  in conjunction with Tcl version 8.6
#  Mar 31, 2024 03:05:40 PM IST  platform: Windows NT
set vTcl(timestamp) ""
if {![info exists vTcl(borrow)]} {
    ::vTcl::MessageBox -title Error -message  "You must open project files from within PAGE."
    exit}


set vTcl(actual_gui_font_dft_desc)  TkDefaultFont
set vTcl(actual_gui_font_dft_name)  TkDefaultFont
set vTcl(actual_gui_font_text_desc)  TkTextFont
set vTcl(actual_gui_font_text_name)  TkTextFont
set vTcl(actual_gui_font_fixed_desc)  TkFixedFont
set vTcl(actual_gui_font_fixed_name)  TkFixedFont
set vTcl(actual_gui_font_menu_desc)  TkMenuFont
set vTcl(actual_gui_font_menu_name)  TkMenuFont
set vTcl(actual_gui_font_tooltip_desc)  TkDefaultFont
set vTcl(actual_gui_font_tooltip_name)  TkDefaultFont
set vTcl(actual_gui_font_treeview_desc)  TkDefaultFont
set vTcl(actual_gui_font_treeview_name)  TkDefaultFont
########################################### 
set vTcl(actual_gui_bg) #d9d9d9
set vTcl(actual_gui_fg) #000000
set vTcl(actual_gui_analog) #ececec
set vTcl(actual_gui_menu_analog) #ececec
set vTcl(actual_gui_menu_bg) #d9d9d9
set vTcl(actual_gui_menu_fg) #000000
set vTcl(complement_color) gray40
set vTcl(analog_color_p) #c3c3c3
set vTcl(analog_color_m) beige
set vTcl(tabfg1) black
set vTcl(tabfg2) white
set vTcl(actual_gui_menu_active_bg)  #ececec
set vTcl(actual_gui_menu_active_fg)  #000000
########################################### 
set vTcl(pr,autoalias) 1
set vTcl(pr,relative_placement) 1
set vTcl(mode) Relative
set vTcl(project_theme) page-dark



proc vTclWindow.top1 {base} {
    global vTcl
    if {$base == ""} {
        set base .top1
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base
    set target $base
    ###################
    # CREATING WIDGETS
    ###################
    vTcl::widgets::core::toplevel::createCmd $top -class Toplevel \
        -menu $top.m73 -background black -highlightbackground black -highlightcolor white 
    wm focusmodel $top passive
    wm geometry $top 601x450+287+156
    update
    # set in toplevel.wgt.
    global vTcl
    global img_list
    set vTcl(save,dflt,origin) 0
    wm maxsize $top 1284 701
    wm minsize $top 120 1
    wm overrideredirect $top 0
    wm resizable $top 1 1
    wm deiconify $top
    set toptitle "Gene Predictor"
    wm title $top $toptitle
    namespace eval ::widgets::${top}::ClassOption {}
    set ::widgets::${top}::ClassOption(-toptitle) $toptitle
    vTcl:DefineAlias "$top" "Toplevel1" vTcl:Toplevel:WidgetProc "" 1
    set vTcl(real_top) {}
    vTcl::widgets::core::popup::createCmd "$top.pop55" \
        -activebackground beige -activeborderwidth 1 -activeforeground black \
        -background black -borderwidth 1 -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground white -tearoff 1 
    global vTcl
    set val vTcl($top.pop55,-proc)
    set vTcl($top.pop55,-proc) popup1
    namespace eval ::widgets::$top.pop55 {}
    set ::widgets::$top.pop55::ClassOption(-proc) popup1
    set ::widgets::$top.pop55::options(-proc) popup1
    set ::widgets::$top.pop55::save(-proc) 1
    vTcl:DefineAlias "$top.pop55" "Popupmenu1" vTcl:WidgetProc "" 1
    ttk::button "$top.tBu58" \
        -text "Choose file" -compound left -cursor fleur 
    vTcl:DefineAlias "$top.tBu58" "TButton2" vTcl:WidgetProc "Toplevel1" 1
    bind $top.tBu58 <<SetBalloon>> {
        set ::vTcl::balloon::%W {Click to select a file from computer}
    }
    label "$top.lab49" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background #000080 -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground white \
        -highlightbackground black -highlightcolor white \
        -text "Upload protein interactions file (.tsv)" 
    vTcl:DefineAlias "$top.lab49" "Label2" vTcl:WidgetProc "Toplevel1" 1
    label "$top.lab62" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background #000080 -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground white \
        -highlightbackground black -highlightcolor white \
        -text "Upload seed protein file (.txt)" 
    vTcl:DefineAlias "$top.lab62" "Label2_1_1" vTcl:WidgetProc "Toplevel1" 1
    radiobutton "$top.rad65" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background black -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9 -weight bold" -foreground white \
        -highlightbackground black -highlightcolor white -justify left \
        -text "Hishigaki Algorithm" -variable "selectedButton" 
    vTcl:DefineAlias "$top.rad65" "Radiobutton1_1" vTcl:WidgetProc "Toplevel1" 1
    radiobutton "$top.rad59" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background black -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9 -weight bold" -foreground white \
        -highlightbackground black -highlightcolor white -justify left \
        -text "Majority Voting Algorithm" -variable "selectedButton" 
    vTcl:DefineAlias "$top.rad59" "Radiobutton1" vTcl:WidgetProc "Toplevel1" 1
    ttk::button "$top.tBu70" \
        -text "Predict candidate genes" -compound left 
    vTcl:DefineAlias "$top.tBu70" "TButton3" vTcl:WidgetProc "Toplevel1" 1
    ttk::button "$top.tBu64" \
        -text "Choose file" -compound left -cursor fleur 
    vTcl:DefineAlias "$top.tBu64" "TButton2_1" vTcl:WidgetProc "Toplevel1" 1
    bind $top.tBu64 <<SetBalloon>> {
        set ::vTcl::balloon::%W {Click to select a file from computer}
    }
    radiobutton "$top.rad68" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background black -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9 -weight bold" -foreground white \
        -highlightbackground black -highlightcolor white -justify left \
        -text "Single function" -variable "selectedButton" 
    vTcl:DefineAlias "$top.rad68" "Radiobutton1_2" vTcl:WidgetProc "Toplevel1" 1
    menu "$top.m73" \
        -activebackground #d9d9d9 -activeforeground black \
        -font "-family {Segoe UI} -size 9" -tearoff 0 
    radiobutton "$top.rad69" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background black -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9 -weight bold" -foreground white \
        -highlightbackground black -highlightcolor white -justify left \
        -text "Multiple functions" -variable "selectedButton" 
    vTcl:DefineAlias "$top.rad69" "Radiobutton1_2_1" vTcl:WidgetProc "Toplevel1" 1
    entry "$top.ent63" \
        -background #c0c0c0 -disabledforeground #3f3f3f \
        -font "-family {Courier New} -size 10" -foreground white \
        -highlightbackground black -highlightcolor white \
        -insertbackground white -selectbackground #d9d9d9 \
        -selectforeground black -width 174 
    vTcl:DefineAlias "$top.ent63" "Entry1_1" vTcl:WidgetProc "Toplevel1" 1
    label "$top.lab74" \
        -activebackground #d9d9d9 -activeforeground #ffffff -anchor w \
        -background black -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9 -slant italic" -foreground #ffffff \
        -highlightbackground black -highlightcolor white \
        -text "Bhagya Wijeratne - 15219" 
    vTcl:DefineAlias "$top.lab74" "Label3" vTcl:WidgetProc "Toplevel1" 1
    label "$top.lab67" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background #000080 -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground white \
        -highlightbackground black -highlightcolor white \
        -text "Single function /  Multiple function" 
    vTcl:DefineAlias "$top.lab67" "Label2_1_2" vTcl:WidgetProc "Toplevel1" 1
    label "$top.lab61" \
        -activebackground #d9d9d9 -activeforeground black -anchor w \
        -background #000080 -compound left -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground #ffffff \
        -highlightbackground black -highlightcolor white \
        -text "Choose Algorithm" 
    vTcl:DefineAlias "$top.lab61" "Label2_1" vTcl:WidgetProc "Toplevel1" 1
    entry "$top.ent57" \
        -background #c0c0c0 -disabledforeground #3f3f3f \
        -font "-family {Courier New} -size 10" -foreground white \
        -highlightbackground black -highlightcolor white \
        -insertbackground white -selectbackground #d9d9d9 \
        -selectforeground black -width 174 
    vTcl:DefineAlias "$top.ent57" "Entry1" vTcl:WidgetProc "Toplevel1" 1
    label "$top.lab47" \
        -activebackground #d9d9d9 -activeforeground #80ffff -anchor w \
        -background #004080 -compound left -disabledforeground #3f3f3f \
        -font "-family {Trebuchet MS} -size 17 -weight bold" \
        -foreground #ffffff -highlightbackground black -highlightcolor white \
        -text "  Network based candidate gene predictor " 
    vTcl:DefineAlias "$top.lab47" "Label1" vTcl:WidgetProc "Toplevel1" 1
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.tBu58 \
        -in $top -x 0 -relx 0.767 -y 0 -rely 0.287 -width 85 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.lab49 \
        -in $top -x 0 -relx 0.117 -y 0 -rely 0.289 -width 0 -relwidth 0.323 \
        -height 0 -relheight 0.047 -anchor nw -bordermode ignore 
    place $top.lab62 \
        -in $top -x 0 -relx 0.117 -y 0 -rely 0.42 -width 0 -relwidth 0.323 \
        -height 0 -relheight 0.047 -anchor nw -bordermode ignore 
    place $top.rad65 \
        -in $top -x 0 -relx 0.518 -y 0 -rely 0.582 -width 0 -relwidth 0.23 \
        -height 0 -relheight 0.056 -anchor nw -bordermode ignore 
    place $top.rad59 \
        -in $top -x 0 -relx 0.518 -y 0 -rely 0.538 -width 0 -relwidth 0.297 \
        -height 0 -relheight 0.056 -anchor nw -bordermode ignore 
    place $top.tBu70 \
        -in $top -x 0 -relx 0.45 -y 0 -rely 0.822 -width 165 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.tBu64 \
        -in $top -x 0 -relx 0.767 -y 0 -rely 0.418 -width 85 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.rad68 \
        -in $top -x 0 -relx 0.518 -y 0 -rely 0.664 -width 0 -relwidth 0.297 \
        -height 0 -relheight 0.056 -anchor nw -bordermode ignore 
    place $top.rad69 \
        -in $top -x 0 -relx 0.518 -y 0 -rely 0.709 -width 0 -relwidth 0.297 \
        -height 0 -relheight 0.056 -anchor nw -bordermode ignore 
    place $top.ent63 \
        -in $top -x 0 -relx 0.467 -y 0 -rely 0.422 -width 174 -relwidth 0 \
        -height 20 -relheight 0 -anchor nw -bordermode ignore 
    place $top.lab74 \
        -in $top -x 0 -relx 0.75 -y 0 -rely 0.933 -width 0 -relwidth 0.34 \
        -height 0 -relheight 0.069 -anchor nw -bordermode ignore 
    place $top.lab67 \
        -in $top -x 0 -relx 0.115 -y 0 -rely 0.691 -width 0 -relwidth 0.323 \
        -height 0 -relheight 0.047 -anchor nw -bordermode ignore 
    place $top.lab61 \
        -in $top -x 0 -relx 0.117 -y 0 -rely 0.558 -width 0 -relwidth 0.323 \
        -height 0 -relheight 0.047 -anchor nw -bordermode ignore 
    place $top.ent57 \
        -in $top -x 0 -relx 0.467 -y 0 -rely 0.289 -width 174 -relwidth 0 \
        -height 20 -relheight 0 -anchor nw -bordermode ignore 
    place $top.lab47 \
        -in $top -x 0 -relx 0.1 -y 0 -rely 0.111 -width 0 -relwidth 0.79 \
        -height 0 -relheight 0.113 -anchor nw -bordermode ignore 

    vTcl:FireEvent $base <<Ready>>
}

proc vTclWindow.top2 {base} {
    global vTcl
    if {$base == ""} {
        set base .top2
    }
    if {[winfo exists $base]} {
        wm deiconify $base; return
    }
    set top $base
    set target $base
    ###################
    # CREATING WIDGETS
    ###################
    vTcl::widgets::core::toplevel::createCmd $top -class Toplevel \
        -menu $top.m87 -background black -highlightbackground black -highlightcolor white 
    wm focusmodel $top passive
    wm geometry $top 600x451+293+124
    update
    # set in toplevel.wgt.
    global vTcl
    global img_list
    set vTcl(save,dflt,origin) 0
    wm maxsize $top 1284 701
    wm minsize $top 120 1
    wm overrideredirect $top 0
    wm resizable $top 1 1
    wm iconify $top
    set toptitle "Results Page"
    wm title $top $toptitle
    namespace eval ::widgets::${top}::ClassOption {}
    set ::widgets::${top}::ClassOption(-toptitle) $toptitle
    vTcl:DefineAlias "$top" "Toplevel2" vTcl:Toplevel:WidgetProc "" 1
    set vTcl(real_top) {}
    vTcl::widgets::core::popup::createCmd "$top.pop75" \
        -activebackground beige -activeborderwidth 1 -activeforeground black \
        -background black -borderwidth 1 -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground white -tearoff 1 
    global vTcl
    set val vTcl($top.pop75,-proc)
    set vTcl($top.pop75,-proc) popup2
    namespace eval ::widgets::$top.pop75 {}
    set ::widgets::$top.pop75::ClassOption(-proc) popup2
    set ::widgets::$top.pop75::options(-proc) popup2
    set ::widgets::$top.pop75::save(-proc) 1
    vTcl:DefineAlias "$top.pop75" "Popupmenu2" vTcl:WidgetProc "" 1
    menu "$top.m87" \
        -activebackground #d9d9d9 -activeforeground black \
        -font "-family {Segoe UI} -size 9" -tearoff 0 
    ttk::label "$top.tLa84" \
        -font "-family {Segoe UI} -size 12 -weight bold" -relief flat \
        -anchor w -justify left -text "Single/ Multiple functions:" \
        -compound left 
    vTcl:DefineAlias "$top.tLa84" "TLabel1_1" vTcl:WidgetProc "Toplevel2" 1
    message "$top.mes85" \
        -background #ffffff -font "-family {Segoe UI} -size 9" \
        -foreground white -highlightbackground black -highlightcolor white \
        -padx 1 -pady 1 -width 140 
    vTcl:DefineAlias "$top.mes85" "Message1_1" vTcl:WidgetProc "Toplevel2" 1
    message "$top.mes96" \
        -background #ffffff -font "-family {Segoe UI} -size 9" \
        -foreground white -highlightbackground black -highlightcolor white \
        -padx 1 -pady 1 -width 60 
    vTcl:DefineAlias "$top.mes96" "Message1_1_2" vTcl:WidgetProc "Toplevel2" 1
    message "$top.mes98" \
        -background #ffffff -font "-family {Segoe UI} -size 9" \
        -foreground white -highlightbackground black -highlightcolor white \
        -padx 1 -pady 1 -width 60 
    vTcl:DefineAlias "$top.mes98" "Message1_1_2_1" vTcl:WidgetProc "Toplevel2" 1
    button "$top.but102" \
        -activebackground #d9d9d9 -activeforeground black -background #004080 \
        -command "Show candidate genes" -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground #ffffff \
        -highlightbackground black -highlightcolor white \
        -text "Predict functions for unknown proteins" 
    vTcl:DefineAlias "$top.but102" "Button1_1" vTcl:WidgetProc "Toplevel2" 1
    button "$top.but101" \
        -activebackground #d9d9d9 -activeforeground black -background #004080 \
        -command "Show candidate genes" -disabledforeground #3f3f3f \
        -font "-family {Segoe UI} -size 9" -foreground #ffffff \
        -highlightbackground black -highlightcolor white \
        -text "Show candidate genes" 
    vTcl:DefineAlias "$top.but101" "Button1" vTcl:WidgetProc "Toplevel2" 1
    ttk::separator "$top.tSe106"
    vTcl:DefineAlias "$top.tSe106" "TSeparator1" vTcl:WidgetProc "Toplevel2" 1
    ttk::label "$top.tLa88" \
        -font "-family {Segoe UI} -size 12 -weight bold" -relief flat \
        -anchor w -justify left -text "Algorithm:" -compound left 
    vTcl:DefineAlias "$top.tLa88" "TLabel1_1_1" vTcl:WidgetProc "Toplevel2" 1
    message "$top.mes89" \
        -background #ffffff -font "-family {Segoe UI} -size 9" \
        -foreground white -highlightbackground black -highlightcolor white \
        -padx 1 -pady 1 -width 260 
    vTcl:DefineAlias "$top.mes89" "Message1_1_1" vTcl:WidgetProc "Toplevel2" 1
    vTcl::widgets::ttk::scrolledtext::CreateCmd "$top.scr105" \
        -borderwidth 2 -relief groove -background black -height 75 \
        -highlightbackground black -highlightcolor white -width 125 
    vTcl:DefineAlias "$top.scr105" "Scrolledtext1" vTcl:WidgetProc "Toplevel2" 1

    $top.scr105.01 configure -background white \
        -font TkTextFont \
        -foreground white \
        -height 3 \
        -highlightbackground black \
        -highlightcolor white \
        -insertbackground white \
        -insertborderwidth 3 \
        -selectbackground #d9d9d9 \
        -selectforeground black \
        -width 10 \
        -wrap none
    button "$top.but103" \
        -activebackground #d9d9d9 -activeforeground black -background #80ffff \
        -disabledforeground #3f3f3f -font "-family {Segoe UI} -size 9" \
        -foreground #004080 -highlightbackground black -highlightcolor white \
        -text "Main Menu" 
    vTcl:DefineAlias "$top.but103" "Button2" vTcl:WidgetProc "Toplevel2" 1
    ttk::label "$top.tLa93" \
        -font "-family {Segoe UI} -size 12 -weight bold" -relief flat \
        -anchor w -justify left -text "No. of proteins in the network:" \
        -compound left 
    vTcl:DefineAlias "$top.tLa93" "TLabel1_1_2" vTcl:WidgetProc "Toplevel2" 1
    ttk::label "$top.tLa97" \
        -font "-family {Segoe UI} -size 12 -weight bold" -relief flat \
        -anchor w -justify left -text "No. of known proteins in the network:" \
        -compound left 
    vTcl:DefineAlias "$top.tLa97" "TLabel1_1_2_1" vTcl:WidgetProc "Toplevel2" 1
    ttk::label "$top.tLa99" \
        -font "-family {Segoe UI} -size 12 -weight bold" -relief flat \
        -anchor w -justify left -text "No. of interactions in the network:" \
        -compound left 
    vTcl:DefineAlias "$top.tLa99" "TLabel1_1_2_1_1" vTcl:WidgetProc "Toplevel2" 1
    message "$top.mes100" \
        -background #ffffff -font "-family {Segoe UI} -size 9" \
        -foreground white -highlightbackground black -highlightcolor white \
        -padx 1 -pady 1 -width 60 
    vTcl:DefineAlias "$top.mes100" "Message1_1_2_1_1" vTcl:WidgetProc "Toplevel2" 1
    ###################
    # SETTING GEOMETRY
    ###################
    place $top.tLa84 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.111 -width 0 -relwidth 0.392 \
        -height 0 -relheight 0.087 -anchor nw -bordermode ignore 
    place $top.mes85 \
        -in $top -x 0 -relx 0.467 -y 0 -rely 0.133 -width 0 -relwidth 0.233 \
        -height 0 -relheight 0.042 -anchor nw -bordermode ignore 
    place $top.mes96 \
        -in $top -x 0 -relx 0.517 -y 0 -rely 0.289 -width 0 -relwidth 0.1 \
        -height 0 -relheight 0.042 -anchor nw -bordermode ignore 
    place $top.mes98 \
        -in $top -x 0 -relx 0.6 -y 0 -rely 0.378 -width 0 -relwidth 0.1 \
        -height 0 -relheight 0.042 -anchor nw -bordermode ignore 
    place $top.but102 \
        -in $top -x 0 -relx 0.55 -y 0 -rely 0.556 -width 217 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.but101 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.556 -width 187 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.tSe106 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.244 -width 0 -relwidth 0.867 \
        -height 0 -relheight 0.004 -anchor nw -bordermode ignore 
    place $top.tLa88 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.044 -width 0 -relwidth 0.175 \
        -height 0 -relheight 0.087 -anchor nw -bordermode ignore 
    place $top.mes89 \
        -in $top -x 0 -relx 0.267 -y 0 -rely 0.067 -width 0 -relwidth 0.433 \
        -height 0 -relheight 0.042 -anchor nw -bordermode ignore 
    place $top.scr105 \
        -in $top -x 0 -relx 0.1 -y 0 -rely 0.644 -width 0 -relwidth 0.837 \
        -height 0 -relheight 0.296 -anchor nw -bordermode ignore 
    place $top.but103 \
        -in $top -x 0 -relx 0.8 -y 0 -rely 0.067 -width 87 -relwidth 0 \
        -height 26 -relheight 0 -anchor nw -bordermode ignore 
    place $top.tLa93 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.267 -width 0 -relwidth 0.408 \
        -height 0 -relheight 0.087 -anchor nw -bordermode ignore 
    place $top.tLa97 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.356 -width 0 -relwidth 0.525 \
        -height 0 -relheight 0.087 -anchor nw -bordermode ignore 
    place $top.tLa99 \
        -in $top -x 0 -relx 0.083 -y 0 -rely 0.443 -width 0 -relwidth 0.458 \
        -height 0 -relheight 0.086 -anchor nw -bordermode ignore 
    place $top.mes100 \
        -in $top -x 0 -relx 0.55 -y 0 -rely 0.466 -width 0 -relwidth 0.1 \
        -height 0 -relheight 0.042 -anchor nw -bordermode ignore 

    vTcl:FireEvent $base <<Ready>>
}

proc 36 {args} {return 1}


Window show .
set btop1 ""
if {$vTcl(borrow)} {
    set btop1 .bor[expr int([expr rand() * 100])]
    while {[lsearch $btop1 $vTcl(tops)] != -1} {
        set btop1 .bor[expr int([expr rand() * 100])]
    }
}
set vTcl(btop) $btop1
Window show .top1 $btop1
if {$vTcl(borrow)} {
    $btop1 configure -background plum
}
set btop2 ""
if {$vTcl(borrow)} {
    set btop2 .bor[expr int([expr rand() * 100])]
    while {[lsearch $btop2 $vTcl(tops)] != -1} {
        set btop2 .bor[expr int([expr rand() * 100])]
    }
}
set vTcl(btop) $btop2
Window show .top2 $btop2
if {$vTcl(borrow)} {
    $btop2 configure -background plum
}

