### elementary steps

# 
try:
    from .groups import fg
except:
    from groups import fg

def build_step(template, groups, **kwargs):
    # Combine groups dict and extra kwargs into one
    all_vars = {**groups, **kwargs}
    return template.format(**all_vars)

def remove_charge(s, charge):
    if charge > 0:
        s = s.replace(f";+{charge}","")
    elif charge == 0:
        try:
            s.replace(f";+0","")
        except:
            return s
    else:
        s = s.replace(f";{charge}","")
    return s

def add_charge(s, charge):
    if charge > 0:
        s += f";+{charge}"
    elif charge == 0:
        s += f";+0"
    else:
        s += f";{charge}"
    return s

def get_charge(s):
    parts = s.split("+")
    if len(parts) >= 2:
        return int(parts[-1][0])
    else:
        parts = s.split("-")
        if len(parts) >= 2:
            try:
                return -int(parts[-1][0])
            except Exception as e:
                # print(e)
                return 0
        else:
            return 0

def change_charge(s, x):
    charge = get_charge(s)
    s = remove_charge(s, charge)
    return add_charge(s, charge + x)

def remove_connection(s):
    parts = s.split("D")
    if len(parts) >= 2:
        h = int(parts[1][0])
        if h > 0:
            h = h - 1
        else:
            raise Exception("No connections to remove")
        parts[1] = f"{h}" + parts[1][1:]
    return "D".join(parts)

def add_connection(s):
    parts = s.split("D")
    if len(parts) >= 2:
        h = int(parts[1][0])
        h += 1
        parts[1] = f"{h}" + parts[1][1:]
    return "D".join(parts)

def deprotonate(s):
    parts = s.split("H")
    if len(parts) >= 2:
        h = int(parts[1][0])
        if h > 0:
            h = h - 1
        else:
            raise Exception("No beta hydrogens to eliminate")
        parts[1] = f"{h}" + parts[1][1:]
    return "H".join(parts)

def protonate(s):
    parts = s.split("H")
    if len(parts) >= 2:
        h = int(parts[1][0])
        h = h + 1
        parts[1] = f"{h}" + parts[1][1:]
    else:
        parts = s.split(":")
        parts[0] += ";H1"
    return "H".join(parts)

steps = {
     
    "Ligand association" : {
        'First': lambda M, lig = "P", c = 0, r = "": 
            build_step("[{lig}{r}:1].[{M};D{c}:2]>>[{lig}:1]-[{M};D{c1}:2]",
                fg,
                M = M,
                lig = lig,
                c = c,
                c1 = c + 1,
                r = r,
            ),
        'Second': lambda M, lig = "P", lig_1 = "P", c = 1, r = "": 
            build_step("[{lig}{r}:1].[{M};D{c}:2]-[{lig_1}:3]>>[{lig}:1]-[{M};D{c1}:2]-[{lig_1}:3]",
                fg,
                M = M, 
                lig = lig,
                lig_1 = lig_1,
                c = c,
                c1 = c + 1,
                r = r,
            ),
        # 'Bis': lambda M, lig = "P", c = 0: 
        #     build_step(f"{fg.get('bis_ligand')(lig, lig, 1, 2)}.[{M};D{c}:3]>>[{M};D{c + 2}:3]([{lig}:1])([{lig}:2])",
        #         fg, 
        #         M = M,
        #         lig = lig,
        #         c = c,
        #     ),
    },

    "Ligand dissociation" : {
        'First': lambda M, lig = "P" : 
            build_step("[{lig}:1]-[{M};D4:2]>>[{lig}:1].[{M};D3:2]",
                fg,
                M = M,
                lig = lig,
            ),
        'Second': lambda M, lig = "P" : 
            build_step("[{lig}:1]-[{M};D3:2]>>[{lig}:1].[{M};D2:2]",
                fg,
                M = M, 
                lig = lig,
            ),
        'Third': lambda M, lig = "P" : 
            build_step("[{lig}:1]-[{M};D2:2]>>[{lig}:1].[{M};D1:2]",
                fg,
                M = M, 
                lig = lig,
            ),
    },

    ### M = metal (Pd, Ni, Cu), l = extra specifications on metal (D0, D2, any exclusions)
    ### l inputs should have leading semicolon
    "Oxidative addition" : (lambda M, l="": 
        build_step("[{C}:1]-[{X}:2].[{M}{l}:3]>>[{C}:1]-[{M}:3]-[{X}:2]", 
        fg, 
        M = M, 
        l = l, 
        )
    ),

    ### M = metal, R1, R2 are combined, l = extra product specifictions (H0, D0, etc.)
    "Reductive elimination" : {
        "Normal" : lambda M, R1, R2, l="": 
            build_step("[{R1}:1]-[{M}:2]-[{R2}:3]>>[{M}{l}:2].[{R1}:1]-[{R2}:3]",
            fg,
            M = M,
            R1 = R1,
            R2 = R2,
            l = l,),
        "Hydride - BA" : lambda M, R, l="": # base-assisted
            build_step("[{M};H1:2]-[{R}:3]>>[{M};H0{l}:2].[{R};H0;-1:3]",
            fg,
            M = M,
            R = R,
            l = l,),
        "Hydride" : lambda M, R, l="": 
            build_step("[{M};H1:2]-[{R}:3]>>[{M};H0{l}:2].[{new_R}:3]",
            fg,
            M = M,
            R = R,
            new_R = protonate(R),
            l = l,),
    },

    ### M1 = metal 1, M2 = metal 2, R1 = group moving from M1 to M2, R2 = group moving from M2 to M1, S = group staying on M1
    ### R1 and R2 must have numbering above 3, R1 and R2 must be bracketed/include full charge information
    "Transmetalation" :  (lambda M1, M2, L1, L2, S="#6":
            build_step("{L1}-[{M1}:1]-[{S}:2].{L2}-[{M2}:3]>>{L1}-[{M2}:3].{L2}-[{M1}:1]-[{S}:2]",
            fg,
            M1 = M1,
            M2 = M2,
            L1 = L1,
            L2 = L2,
            S = S,
            )
    ),

    ### M = metal, A (2) = whatever atom is alpha to M, B (3) is the carbon containing the beta C-H bond for elimination
    "Beta-hydride elimination" : (lambda M, A, B, l = "":
        build_step("[{M}{l}:1]-[{A}:2]-[{B}:3]>>[{M};H1:1].[{A}:2]=[{B_h}:3]",
        fg, 
        M = M,
        A = A,
        B = B,           
        B_h = deprotonate(B),
        l = l,
        )
    ),

    ## beta carbon elimination

    ### M = metal, R = R-group that is migrating, A - B are bonded by either double or triple bond set by bond type
    "Migratory insertion" : {
        
        'Normal': (lambda M, R, A, B, b_type = "=", l = "":
            build_step("[{R}:2]-[{M}{l}:1].[{A}:3]{b_type}[{B}:4]>>[{R}:2]-[{A}:3]{new_b}[{B}:4]-[{M}{l}:1]",
            fg, 
            M = M,
            A = A,
            B = B,
            R = R,
            b_type = b_type,
            new_b = "=" if b_type == "#" else "-",
            l = l,
            )
        ),
        'Hydride' : (lambda M, num_H, A, B, b_type = "=", l = "":
            build_step("[{M};H{num_H}{l}:1].[{A}:3]{b_type}[{B}:4]>>[{new_A}:3]{new_b}[{B}:4]-[{M};H{new_H}{l}:1]",
            fg, 
            M = M,
            num_H = num_H,
            new_H = num_H - 1,
            A = A,
            new_A = protonate(A),
            B = B,
            b_type = b_type,
            new_b = "=" if b_type == "#" else "-",
            l = l,
            )
        ),
    },

    ### hydroxy group exchanged with halide on metal
    "Halide exchange" : {
        # "Water": lambda M,: 
        #     build_step("[{H2O}:1].[{X}:2]-[{M}:3]-[#6:4]>>[{X};H0;-1:2].[{OH}:1]-[{M}:3]-[#6:4]", 
        #     fg,
        #     M = M,),
        "Hydroxide" : lambda M, l = "-[#6:4]": 
            build_step("[{OH_}:1].[{X}:2]-[{M}:3]{l}>>[{X};-1:2].[{OH}:1]-[{M}:3]{l}", 
            fg,
            M = M,
            l = l,),
        "Hydroxide-Metal" : lambda M, l = "-[#6:4]": 
            build_step("[{OH}:1]-[{AM}:5].[{X}:2]-[{M}:3]{l}>>[{AM};+1:5].[{X};-1:2].[{OH}:1]-[{M}:3]{l}", 
            fg,
            M = M,
            l = l,),
        # "Alcohol" : lambda M,: 
        #     build_step("[{not_carbonyl}:1]-[{OH}:5].[{X}:2]-[{M}:3]-[#6:4]>>[{X};H0;-1:2].[{C}:1]-[{O_}:5]-[{M}:3]-[#6:4]", 
        #     fg,
        #     M = M,),
        "Alkoxide" : lambda M, l = "-[#6:4]": 
            build_step("[{not_carbonyl}:1]-[{O_};-1:5].[{X}:2]-[{M}:3]{l}>>[{X};-1:2].[{C}:1]-[{O_};+0:5]-[{M}:3]{l}", 
            fg,
            M = M,
            l = l,),
        "Alkoxide-Metal" : lambda M, l = "-[#6:4]": 
            build_step("[{not_carbonyl}:1]-[{O_}:5]-[{AM}:6].[{X}:2]-[{M}:3]{l}>>[{AM};+1:6].[{X};-1:2].[{C}:1]-[{O_};+0:5]-[{M}:3]{l}", 
            fg,
            M = M,
            l = l,),
        # "Carbonate base" :,

        #### FIX HOW THIS WORKS, ALCOHOL AND WATER SHOULD BE TWO DISCRETE STEPS
    },

    "PdOB formation" : {
        "Pd-hydroxy" : "[#6:3]-[Pd:1]-[#8:2].[#5;+0:4]-[#6;+0:5]>>[#6:3]-[Pd:1]-[#8;+1:2]-[#5;-1:4]-[#6;+0:5]",
        "Borate" : f"[#6:3]-[Pd:1]-[{fg.get('X')}:2].[#5;-1:4]-[#8;+0;!$([#8]-[#6](-[#6])(-[#6])-[#6](-[#6])(-[#6])):5]>>[{fg.get('X')};-1:2].[#6:3]-[Pd:1]-[#8;+1:5]-[#5;-1:4]"},

    "Borate formation" : {
        "Hydroxide" : "[#5;!$([#5]-[#9]):2]-[#6;+0:4].[O;-1;H1:1]>>[#6;+0:4]-[#5;-1:2]-[O;+0;H1:1]",
        "Metal-hydroxide/alkoxide" : f"[#5;!$([#5]-[#9]):2]-[#6;+0:4].[O;$([O]-[C,#1]);!$([O]-[C]=[O]):1]-[{fg.get('AM')}:3]>>[#6;+0:4]-[#5;-1:2]-[O;+0:1].[{fg.get('AM')};+1:3]",
        "Alkoxide" : f"[O;$([O]-[C]);!$([O]-[C]=[O]);-1;H0:1].[#5;!$([#5]-[#9]):2]-[#6;+0:4]>>[#6;+0:4]-[#5;-1:2]-[O;+0:1]",
    },

    "Fluoride leaving" : {
        "First" : "[#5;$([#5](-[#9])(-[#9])(-[#9]));-1:1]-[#9:2]>>[#5;+0:1].[#9;-1:2]",
        "Second" : "[#5;$([#5](-[#9])(-[#9]));-1:1]-[#9:2]>>[#5;+0:1].[#9;-1:2]",
        "Third" : "[#5;$([#5](-[#9]));-1:1]-[#9:3]>>[#5;+0:1].[#9;-1:3]",
    },
    
    "Phosphine reduction" : {
        "X" : lambda M, lig :
            build_step("[{X}:3]-[{M}:1]-[{X}:2].[{lig}:4]>>[{X}:2]-[{lig};+1:4].[{X};-1:3].[{M}:1]",
            fg,
            M = M, 
            lig = lig,
            )
    },
    
    "SN2" : ( lambda nuc, elec, lg, l = "" : 
        build_step("[{elec}{l}:1]-[{lg}:2].[{nuc}:3]>>[{charged_nuc}:3]-[{elec}:1].[{charged_lg}:2]",
        fg,
        nuc = nuc, 
        elec = elec,
        lg = lg,
        l = l,
        charged_nuc = change_charge(nuc, 1),
        charged_lg = change_charge(lg, -1),
        )
    ),

    # 'PtO2 reduction' : '[O:2]=[Pt:1]=[O:3].[H:4]-[H:5].[H:6]-[H:7]>>[H]-[O:2]-[H].[H]-[O:3]-[H].[Pt:1]'
    



}  
