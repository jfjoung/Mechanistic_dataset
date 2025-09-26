## for running template creation

try:
    from .groups import fg
    from .elem_steps import *
except:
    from groups import fg
    from elem_steps import *

import re

# def more_subs(text):
#     def add_two(match):
#         number = int(match.group(1))
#         return f";D{number + 2}"
    
#     return re.sub(r";D(\d+)", add_two, text)

def more_subs(text):
    match = re.search(r";D(\d+)", text)
    if match:
        number = int(match.group(1))
        return f";D{number + 2}"
    return None 


class AttrDict(dict):
    def __getattr__(self, key):
        return self.get(key)

fg = AttrDict(fg)

# M1 = metal 1 (catalyst in entire cycle), M2 = metal that R group is transmetalating from, i.e. B in suzuki, L1 = group transmetalating from M1 to M2, L2 = group transmetalating from M2 to M1,
# RE1 and RE2 two groups that undergo reductive elimination to form product
def core_cat_cycle (stage_nums, M1, M2, L1, L2, RE1, RE2, l1="", l2=""):
    cycle = {}

    if stage_nums[0] is not None:
        cycle[stage_nums[0]] = {
            'Templates': [steps["Oxidative addition"](M1, l1),],
            'pKa': [None],
            'Description': 'Oxidative addition'
        }
    if stage_nums[1] is not None:
        cycle[stage_nums[1]] = {
            'Templates': [steps["Transmetalation"](M1, M2, L1, L2)],
            'pKa': [None],
            'Description': 'Transmetalation'
        }
    if stage_nums[2] is not None:
        cycle[stage_nums[2]] = {
            'Templates': [steps["Reductive elimination"]["Normal"](M1, RE1, RE2, l2)],
            'pKa': [None],
            'Description': 'Reductive elimination'
        }

    return cycle

common_steps = {
    # 'Halide exchange with hydroxide' : lambda offset, M, :
    # {
        
    # }
}

catalyst_preparation = {
    # offset = step num. to start catalyst preparation at (usually 0), M = metal, lig = 1st ligand to detach/attach
    # ligand_2 = 2nd ligand to detach/attach, defaults to lig, c = number of connections on the metal before dis/association
    "Attach ligands" : lambda offset, M, lig, lig_restriction = "" ,c = 0, lig_2 = None, :
        {offset : {
            'Templates': [steps["Ligand association"]["First"](M, lig, c, lig_restriction),],
            'pKa': [None],
            'Description': 'Attach first ligand'
        },
        offset + 1 : {
            'Templates': [steps["Ligand association"]["Second"](M, (lig_2 if lig_2 is not None else lig), lig, c+1, lig_restriction),],
            'pKa': [None],
            'Description': 'Attach second ligand'
        },
        },
    "Remove ligands" : lambda offset, M, lig:
        {offset : {
            'Templates': [steps["Ligand dissociation"]["First"](M, lig),],
            'pKa': [None],
            'Description': 'Remove fourth ligand'
        },
        offset + 1 : {
            'Templates': [steps["Ligand dissociation"]["Second"](M, lig),],
            'pKa': [None],
            'Description': 'Remove third ligand'
        }},
    "Remove three ligands" : lambda offset, M, lig:
        {offset : {
            'Templates': [steps["Ligand dissociation"]["First"](M, lig),],
            'pKa': [None],
            'Description': 'Remove fourth ligand'
        },
        offset + 1 : {
            'Templates': [steps["Ligand dissociation"]["Second"](M, lig),],
            'pKa': [None],
            'Description': 'Remove third ligand'
        },
        offset + 2 : {
            'Templates': [steps["Ligand dissociation"]["Third"](M, lig),],
            'pKa': [None],
            'Description': 'Remove second ligand'
        }
        },
    "PdX2 with phosphines water and base" : {
        'Monophosphine' : 
            lambda offset, M, lig:
                {offset : {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
                }, 
                offset + 1:{
                    'Templates': [steps['Halide exchange']["Hydroxide"]("Pd", ""),
                                  steps['Halide exchange']["Hydroxide-Metal"]("Pd", ""),],
                    'pKa': [None, None],
                    'Description': 'Halide exchange on Pd(II)'
                },
                offset + 2: {
                    'Templates': [f"[{lig}:1]-[{M};D4:2]-[O;H1:3]>>([{M};-1:2].[O;H1:3]-[{lig};+1:1])"],
                    'pKa': [None],
                    'Description': 'Reductive elimination'
                },
                offset + 3 : {
                    'Templates': ["[P;H0;+1:4]-[OH1;+0:5]>>[P;H0;+0:4]=[O;H0;+0:5]"],
                    'pKa': [{'B': -5}],
                    'Description': 'Formation of phosphine oxide'
                },
                offset + 4 : {
                    'Templates': [f"[{M};-1:1]-[{fg.X}:2].[{lig};!$([P]-[Pd]):3]>>[{M};+0:1]-[{lig}:3].[{fg.X};-1:2]"],
                    'pKa': [None],
                    'Description' : 'Halide exchange'
                },
                offset + 5 : {
                    'Templates': [f"[{M};-1:1]-[{fg.X}:2]>>[{M};+0:1].[{fg.X};-1:2]"],
                    'pKa': [None],
                    'Description' : 'Halide leaves'
                }
                },
        'Monophosphine - 2 eq' : 
            lambda offset, M, lig:
                {offset : {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
                }, 
                offset + 1:{
                    'Templates': [steps['Halide exchange']["Hydroxide"]("Pd", "")],
                    'pKa': [None],
                    'Description': 'Halide exchange on Pd(II)'
                },
                offset + 2: {
                    'Templates': [f"[{lig}:1]-[{M};D4:2]-[O;H1:3]>>[{M};-1:2].[O;H1:3]-[{lig};+1:1]"],
                    'pKa': [None],
                    'Description': 'Reductive elimination'
                },
                offset + 3 : {
                    'Templates': ["[P;H0;+1:4]-[OH1;+0:5]>>[P;H0;+0:4]=[O;H0;+0:5]"],
                    'pKa': [{'B': -5}],
                    'Description': 'Formation of phosphine oxide'
                },
                offset + 4 : {
                    'Templates': [f"[{M};-1:1]-[{fg.X}:2]>>[{M};+0:1].[{fg.X};-1:2]"],
                    'pKa': [None],
                    'Description' : 'Halide leaves'
                }
                },
        'Bisphosphine' : 
            lambda offset, M, lig:
                {offset : {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
                }, 
                offset + 1:{
                    'Templates': [steps['Halide exchange']["Hydroxide"]("Pd", "")],
                    'pKa': [None],
                    'Description': 'Halide exchange on Pd(II)'
                },
                offset + 2: {
                    'Templates': [f"[{lig}:4]1-[*:5]-[*:6]-[*:7]-[{lig}:1]-[{M};D4:2]1-[O;H1:3]>>[O;H1:3]-[{lig};+1:1]-[*:5]-[*:6]-[*:7]-[{lig}:4]-[{M};-1:2]"],
                    'pKa': [None],
                    'Description': 'Reductive elimination'
                },
                offset + 3 : {
                    'Templates': ["[P;H0;+1:4]-[OH1;+0:5]>>[P;H0;+0:4]=[O;H0;+0:5]"],
                    'pKa': [{'B': -5}],
                    'Description': 'Formation of phosphine oxide'
                },
                offset + 4 : {
                    'Templates': [f"[{M};-1:1]-[{fg.X}:2].[{lig};!$([P]-[Pd]):3]>>[{M};+0:1]-[{lig}:3].[{fg.X};-1:2]"],
                    'pKa': [None],
                    'Description' : 'Halide exchange'
                }, 
                offset + 5 : {
                    'Templates': [steps["Ligand association"]["Second"](M, lig, lig, 2, "!$([P]-[Pd])"),],
                    'pKa': [None],
                    'Description': 'Attach second part of bisphosphine'
                },
                offset + 6 : {
                    'Templates': [f"[{M};D3:1]-[{lig};$({fg.p_p_oxide}):2]>>[{M};+0:1].[{lig}:2]"],
                    'pKa': [None],
                    'Description' : 'Dissociate oxidized ligand'
                },
                },
        'Bisphosphine - 1 eq' : 
            lambda offset, M, lig:
                {offset : {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
                }, 
                offset + 1:{
                    'Templates': [steps['Halide exchange']["Hydroxide"]("Pd", "")],
                    'pKa': [None],
                    'Description': 'Halide exchange on Pd(II)'
                },
                offset + 2: {
                    'Templates': [f"[{lig}:4]1-[*:5]-[*:6]-[*:7]-[{lig}:1]-[{M};D4:2]1-[O;H1:3]>>[O;H1:3]-[{lig};+1:1]-[*:5]-[*:6]-[*:7]-[{lig}:4]-[{M};-1:2]"],
                    'pKa': [None],
                    'Description': 'Reductive elimination'
                },
                offset + 3 : {
                    'Templates': ["[P;H0;+1:4]-[OH1;+0:5]>>[P;H0;+0:4]=[O;H0;+0:5]"],
                    'pKa': [{'B': -5}],
                    'Description': 'Formation of phosphine oxide'
                },
                offset + 4 : {
                    'Templates': [f"[{M};-1:1]-[{fg.X}:2].[O;$([O]=[P]-[#6]):3]>>[{M};+0:1]-[O:3].[{fg.X};-1:2]"],
                    'pKa': [None],
                    'Description' : 'Halide exchange'
                }, 
                },
    },
    "PdX2 with phosphines and no water" : {
        'Monophosphine' : 
            lambda offset, M, lig:
                {offset : {
                    'Templates': [f"[{fg.X}:1]-[{M};D4:2](-[{lig}:4])-[{fg.X}:3]>>([{lig};+1:4]-[{fg.X}:3].[{M}:2].[{fg.X};-1:1])"],
                    'pKa': [None], 
                    'Description': 'Reduction of Pd(II)'
                }, 
            },
        # 'Bisphosphine' : 
        #     lambda offset, M, lig:
        #             {offset : {
        #                 'Templates': [f"[{fg.X}:1]-[{M}:2](-[{lig}:4])-[{fg.X}:3]>>([{lig};+1:4]-[{fg.X}:3].[{M}:2].[{fg.X};-1:1])"],
        #                 'pKa': [None], 
        #                 'Description': 'Reduction of Pd(II)'
        #             }, 
        #         },
    },
    "PdX2 with amine base" : lambda offset, M = "Pd", B = "N;!$([N]-[C]=[O]);$([N]-[C]-[#1])", l = "" : 
        {offset : {
            'Templates': [steps["SN2"](B + ";+0", M, fg.X,l)],
            'pKa': [None],
            'Description': 'Ligand exchange'
        },
        offset + 1 : {
            'Templates': [steps["Beta-hydride elimination"](M, B + ";+1", "#6;H2"),
                          steps["Beta-hydride elimination"](M, B + ";+1", "#6;H1"),],
            'pKa': [None, None],
            'Description': 'Beta-hydride elimination'
        },
        offset + 2: {
            'Templates': [steps["Reductive elimination"]["Hydride - BA"](M, fg.X)],
            'pKa': [{'B': -6}],
            'Description': 'Reductive elimination',
        },},
    "PdX2 with phosphines and phosphate/carbonate" : lambda offset, lig:
        {offset : {
            'Templates': [f"[{fg.anion}:7]-[Pd;D4:1].[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O])):2]-[{fg.AM}:3]>>[{fg.anion};-1:7].[Pd:1]-[O:2].[{fg.AM};+1:3]",
                          f"[{fg.anion}:7]-[Pd;D4:1].[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O]));-1:2]>>[{fg.anion};-1:7].[Pd:1]-[O;+0:2]"],
            'pKa': [None, None],
            'Description': 'Halide exchange',
        },
        offset + 1 : {
            'Templates': [f"[P:2]-[Pd;D4:1]-[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O])):3]>>([P;+1:2]-[O:3].[Pd;D2;-1:1])"],
            'pKa': [None],
            'Description': 'POAc+ leaving'
        },
        offset + 2 : {
            'Templates': [f"[P;+1:1]-[O:2]-[C,P:3]=[O:4].[O;!$([O]-[P,C]-[O]-[P;+1]);$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O])):5]-[{fg.AM}:6]>>[P;+1:1]-[O:2]-[C,P:3](-[O;-1:4])-[O:5].[{fg.AM};+1:6]",
                          "[P;+1:1]-[O:2]-[C,P:3]=[O:4].[O;!$([O]-[P,C]-[O]-[P;+1]);$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O]));-1:5]>>[P;+1:1]-[O:2]-[C,P:3](-[O;-1:4])-[O;+0:5]"],
            'pKa': [None, None],
            'Description' : 'Carbonate/Phosphate attack'
        },
        offset + 3: {
            'Templates' : ["[P;+1:1]-[O:2]-[C,P:3]-[O;-1:4]>>[P;+0:1]=[O:2].[O;+0:4]=[C,P:3]"],
            'pKa' : [None],
            'Description' : "P=O formation"
        },
        offset + 4:{
            'Templates': [steps["Ligand association"]["Second"]("Pd", lig, lig, 2),],
            'pKa': [None],
            'Description': 'Attach second ligand'
        },
        offset + 5: {
            'Templates' : [f"[Pd;-1:1]-[{fg.anion}:2]>>[Pd;+0:1].[{fg.anion};-1:2]"],
            'pKa' : [None],
            'Description' : 'Halide leaves'
        }
        },
    "Pd-allyl" : lambda offset : 
        {offset : {
            'Templates' : [steps['Transmetalation']("Pd;D3", "Zn,Sn,B,Mg,Pd", f"[{fg.X}:4]", "[#6:5]")],
            'pKa' : [None],
            'Description' : 'Transmetalation'
        },
        offset + 1 : {
            'Templates' : [steps['Reductive elimination']['Normal']("Pd", "#6", "#6")],
            'pKa' : [None],
            'Description' : 'Reductive elimination'
        }
        },
    "Reduction of PdX2" : lambda offset : 
        {offset : {
            'Templates' : [f"[Pd;D4:1](-[{fg.anion}:2])(-[{fg.anion}:3]).[Zn,Sn,Mg;D0:4]>>[Pd;D2:1].[Zn,Sn,Mg;D2:4](-[{fg.anion}:2])(-[{fg.anion}:3])"],
            'pKa' : [None],
            'Description' : 'Reduction of Pd'
        },
        },
    "Buchwald precatalyst" : lambda offset, : 
        {offset : {
            'Templates': [f"[{fg.anion}:7]-[Pd:1]-[#6:2]:[#6:3]-[#6:4]:[#6:5]-[#7;H2:6]>>[{fg.anion};-1:7].[Pd:1]1-[#6:2]:[#6:3]-[#6:4]:[#6:5]-[#7;H1:6]-1",
                          f"[{fg.anion}:7]-[Pd:1]-[#6:2]:[#6:3]-[#6:4]:[#6:5]-[#7;H1:6]>>[{fg.anion};-1:7].[Pd:1]1-[#6:2]:[#6:3]-[#6:4]:[#6:5]-[#7;H0:6]-1"],
            'pKa': [{'B': 3}, {'B': 3}],
            'Description': 'Deprotonation of precatalyst'
        },
        offset + 1 : {
            'Templates': [f"[Pd:1]1-[#6:2]:[#6:3]-[#6:4]:[#6:5]-[#7:6]-1>>[#6:2]1:[#6:3]-[#6:4]:[#6:5]-[#7:6]-1.[Pd:1]"],
            'pKa': [None],
            'Description': 'Reductive elimination to Pd(0)'
        },},
    "Pd(OAc)2 with phosphines" : 
        {'Monophosphine' : lambda offset, lig : 
            {offset : {
                'Templates': [f"[P:2]-[Pd;D4:1](-[O:6]-[C:7]=[O:8])-[O:3]-[C:4]=[O:5]>>[P;+1:2]-[O:5]-[C:4]=[O:3].[Pd;D2;-1:1]-[O:6]-[C:7]=[O:8]"],
                'pKa': [None],
                'Description': 'POAc+ leaving'
            },
            offset + 1: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
            }, 
            offset + 2: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3]=[O:4].[O;H1;-1:5]>>[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;H1;+0:5]"],
                'pKa': [None],
                'Description': 'Hydrolysis to form phosphine oxide'
            },
            offset + 3: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;H1;+0:5]>>[P;+0:1]=[O:2].[C:3](=[O;+0:4])-[O;H1;+0:5]"],
                'pKa': [None],
                'Description': 'P=O leaving'
            },
            offset + 4: {
                'Templates': [f"[Pd;D2;-1:1]-[O:2]-[C:3]=[O:4].[{lig};!$([P]-[Pd]):5]>>[{lig}:5]-[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'Ligand association'
            },
            offset + 5: {
                'Templates': [f"[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]>>[Pd;D2;+0:1].[O;-1:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'OAc- leaves'
            }
            },
        'Bisphosphine' : lambda offset : 
            {offset : {
                'Templates': [f"[P:2]-[Pd;D4:1](-[O:6]-[C:7]=[O:8])-[O:3]-[C:4]=[O:5]>>[P;+1:2]-[O:5]-[C:4]=[O:3].[Pd;D2;-1:1]-[O:6]-[C:7]=[O:8]"],
                'pKa': [None],
                'Description': 'POAc+ leaving'
            },
            offset + 1: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]"],
                    'pKa': [{'B': -2}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide from Water'
            }, 
            offset + 2: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3]=[O:4].[O;H1;-1:5]>>[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;H1;+0:5]"],
                'pKa': [None],
                'Description': 'Hydrolysis to form phosphine oxide'
            },
            offset + 3: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;H1;+0:5]>>[P;+0:1]=[O:2].[C:3](=[O;+0:4])-[O;H1;+0:5]"],
                'pKa': [None],
                'Description': 'P=O leaving'
            },
            offset + 4: {
                'Templates': [f"[Pd;D2;-1:1]-[O:2]-[C:3]=[O:4].[O;$([O]=[P]-[#6]):5]>>[O:5]-[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'Ligand association'
            },
            offset + 5: {
                'Templates': [f"[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]>>[Pd;D2;+0:1].[O;-1:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'OAc- leaves'
            }
            },

        },
     "Pd(OAc)2 without water" : 
        {'Monophosphine' : lambda offset, lig : 
            {offset : {
                'Templates': [f"[P:2]-[Pd;D4:1](-[O:6]-[C:7]=[O:8])-[O:3]-[C:4]=[O:5]>>[P;+1:2]-[O:5]-[C:4]=[O:3].[Pd;D2;-1:1]-[O:6]-[C:7]=[O:8]"],
                'pKa': [None],
                'Description': 'POAc+ leaving'
            },
            offset + 1: {
                'Templates': [f"[Pd;D2;-1:1]-[O:2]-[C:3]=[O:4].[{lig};!$([P]-[Pd]):5]>>[{lig}:5]-[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'Ligand association'
            },
            offset + 2: {
                'Templates': [f"[Pd;D3;-1:1]-[O:2]-[C:3]=[O:4]>>[Pd;D2;+0:1].[O;-1:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'OAc- leaves'
            },
            offset + 3: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3]=[O:4].[O;-1:5]-[C:7]=[O:8]>>[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;+0:5]-[C:7]=[O:8]"],
                'pKa': [None],
                'Description': 'OAc- attack'
            },
            offset + 4: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;+0:5]>>[P;+0:1]=[O:2].[C:3](=[O;+0:4])-[O;+0:5]"],
                'pKa': [None],
                'Description': 'P=O leaving'
            },
            },
        'Bisphosphine' : lambda offset, lig : 
            {offset : {
                'Templates': [f"[P:2]-[Pd;D4:1](-[O:6]-[C:7]=[O:8])-[O:3]-[C:4]=[O:5]>>([P;+1:2]-[O:5]-[C:4]=[O:3].[Pd;D2;-1:1]-[O:6]-[C:7]=[O:8])"],
                'pKa': [None],
                'Description': 'POAc+ leaving'
            },
            offset + 1: {
                'Templates': [f"[Pd;D2;-1:1]-[O:2]-[C:3]=[O:4]>>[Pd;D1;+0:1].[O;-1:2]-[C:3]=[O:4]"],
                'pKa': [None],
                'Description': 'OAc- leaves'
            },
            offset + 2: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3]=[O:4].[O;-1:5]-[C:7]=[O:8]>>[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;+0:5]-[C:7]=[O:8]"],
                'pKa': [None],
                'Description': 'OAc- attack'
            },
            offset + 3: {
                'Templates': [f"[P;+1:1]-[O:2]-[C:3](-[O;-1:4])-[O;+0:5]>>[P;+0:1]=[O:2].[C:3](=[O;+0:4])-[O;+0:5]"],
                'pKa': [None],
                'Description': 'P=O leaving'
            },
            },

        }
}

def general_cu_coupling (cycle, offset=0, global_reagent=[], global_exclude=[], modifier="", partner=""):
    return {
        f'Cu(I)-X ligandless {modifier}':{
            'Reagent': [] + global_reagent,
            'Exclude_reagent': [f'[{fg.n_n_ligand_check}]',f'[{fg.n_o_ligand}]','[O]=[C]-[C]-[C]=[O]','[Pd]'],
            'Stages': cycle(offset, partner , 1)
        },
        f'Cu(I)-X with ligand {modifier}':{
            'Reagent': [f'[{fg.n_n_ligand_check}]'] + global_reagent,
            'Exclude_reagent': ['[Pd]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Cu", f"#7;{fg.n_n_ligand}", c = 1),
                **cycle(offset + 2, partner , 3)}
        },
        f'Cu(I)-X with n,o ligand {modifier}':{
            'Reagent': [f'[{fg.n_o_ligand}]'] + global_reagent,
            'Exclude_reagent': ['[Pd]'],
            'Stages': {
                offset : {
                    'Templates': ["[O;H1:1]-[C:2](=[O:3])-[C:4]-[N:5]>>[O;-1;H0:1]-[C:2](=[O:3])-[C:4]-[N:5]",
                                  "[O;H1:1]-[c:2]:[c:4]:[n:5]>>[O;H0;-1:1]-[c:2]:[c:4]:[n:5]",
                                  "[O;H1:1]-[C:2](=[O:3])-[c:4]:[n:5]>>[O;H0;-1:1]-[C:2](=[O:3])-[c:4]:[n:5]"
                                  ],
                    'pKa': [{'B' : 3}, {'B' : 3}, {'B' : 3}],
                    'Description': 'Deprotonate'
                },
                offset + 1 : {
                    'Templates': ["[O;$([O;-1]-[C](=[O])-[C]-[N]),$([O;-1]-[c]:[c]:[n]),$([O;-1]-[C](=[O])-[c]:[n]),$([O;-1]-[C](=[O])-[c]:[n]):1].[Cu;D1:2]>>[Cu;-1:2]-[O;+0:1]"],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                offset + 2 : {
                    'Templates': ["[Cu;-1:1]-[#8:2].[#7;$([N]-[C]-[C](=[O])-[O]),$([n]:[c]:[c]-[O]),$([n]:[c]-[C](=[O])-[O]),$([n]:[c]-[C](=[O])-[O]):3]>>[#7:3]-[Cu:1]-[#8:2]"],
                    'pKa': [None],
                    'Description': 'Attach second ligand'
                },
                **cycle(offset + 3, partner, 3)}
        },
        f'Cu(I)-X with diketone {modifier}':{
            'Reagent': ['[O]=[C]-[C]-[C]=[O]'] + global_reagent,
            'Exclude_reagent': ['[Pd]'],
            'Stages': {
                offset : {
                    'Templates': ["[O:1]=[C:2]-[C;H1:3]-[C:4]=[O:5]>>[O;-1:1]-[C:2]=[C;H0:3]-[C:4]=[O:5]",
                                  "[O:1]=[C:2]-[C;H2:3]-[C:4]=[O:5]>>[O;-1:1]-[C:2]=[C;H1:3]-[C:4]=[O:5]"
                                  ],
                    'pKa': [{'B' : 3},{'B' : 3}],
                    'Description': 'Deprotonate'
                },
                offset + 1 : {
                    'Templates': ["[O;-1;$([O]-[C]=[C]-[C]=[O]):1].[Cu;D1:2]>>[Cu;-1:2]-[O;+0:1]"],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                offset + 2 : {
                    'Templates': ["[Cu;-1:1]-[#8:2].[O;$([O]=[C]-[C]=[C]-[O]-[Cu]):3]>>[#8:3]-[Cu:1]-[#8:2]"],
                    'pKa': [None],
                    'Description': 'Attach second ligand'
                },
                **cycle(offset + 3, partner, 3)}
        },
    }

def general_pd_coupling (cycle, offset=0, global_reagent=[], global_exclude=[], modifier=""):
    return {
        f'Pd(0) ligandless {modifier}': {
            'Reagent': ['[Pd;D0;+0]'] + global_reagent,
            'Exclude_reagent': [f'[{fg.p_ligand}]'] + global_exclude,
            'Stages': cycle(offset, ";D0"),
        },
        f'Pd(0) with ligand - bis {modifier}': {
            'Reagent': ['[Pd;D0;+0]', f'[{fg.p_ligand}]'] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])"),
                **cycle(offset + 2, ";D2;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'Pd(0) with ligand - mono {modifier}': {
            'Reagent': ['[Pd;D0;+0]', f'[{fg.p_ligand}]'] + global_reagent,
            'Exclude_reagent': [ fg.bisphosphine] + global_exclude,
            'Stages': {
                offset : {
                    'Templates': [steps["Ligand association"]["First"]("Pd", fg.p_ligand, 0, "!$([P]-[Pd])"),],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                **cycle(offset + 1, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'Pd with two ligands {modifier}': {
            'Reagent': [f'[Pd;D2](-[P])-[P]'] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **cycle(offset, ";D2;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'Pd with three or four ligands - bis {modifier}' : {
            'Reagent': [f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **catalyst_preparation['Remove ligands'](offset, "Pd", fg.p_ligand),
                **cycle(offset + 2, ";D2;!$([Pd]-[Cl,Br,I,O])")},
        },
        f'Pd(0) with two, three, or four ligands {modifier}' : {
            'Reagent': [f'[Pd;D2,D3,D4;$([Pd](-[*])-[*]),$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'] + global_reagent,
            'Exclude_reagent': [fg.bisphosphine] + global_exclude,
            'Stages': {
                **catalyst_preparation['Remove three ligands'](offset, "Pd", fg.p_ligand),
                **cycle(offset + 3, ";D1;!$([Pd]-[Cl,Br,I,O])")
            },
        },
        f'Pd Buchwald precatalyst {modifier}' : {
            'Reagent' : [f'[{fg.anion}]-[Pd](-[P;H0])-[#6]:[#6]-[#6]:[#6]-[#7;H2,H1]'] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **catalyst_preparation['Buchwald precatalyst'](offset),
                **cycle(offset + 2, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },

        f'PdX2, L2 with amine {modifier}': {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', fg.amine] + global_reagent,
            'Exclude_reagent': ['[O;H2]'] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with amine base'](offset + 2, "Pd;D4"),
                **cycle(offset + 5, ";D2;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2, monophosphine with water {modifier}': {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', fg.oh_or] + global_reagent,
            'Exclude_reagent': [fg.bisphosphine] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine'](offset + 2, "Pd", fg.p_ligand),
                **cycle(offset + 8, ";D2,D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2, bisphosphine with water {modifier}': {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', fg.oh_or, fg.bisphosphine] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Bisphosphine'](offset + 2, "Pd", fg.p_ligand),
                **cycle(offset + 9, ";D2,D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2, L2 with phosphate/carbonate {modifier}': {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', '[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O]))]'] + global_reagent,
            'Exclude_reagent': ['[O;H2]'] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines and phosphate/carbonate'](offset + 2, fg.p_ligand),
                **cycle(offset + 8, ";D2,D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2, L2 {modifier}' : {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]'] + global_reagent,
            'Exclude_reagent': ['[O;H2]', '[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O]))]', fg.amine, fg.bisphosphine] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines and no water']['Monophosphine'](offset + 2, "Pd", fg.p_ligand),
                **cycle(offset + 3, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2, L2 reduction {modifier}' : {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', '[Zn,Sn,Mg;D0]'] + global_reagent,
            'Exclude_reagent': ['[O;H2]', '[O;$([O]-[P](=[O])(-[O])(-[O])),$([O]-[C](=[O])(-[O]))]', fg.amine] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Reduction of PdX2'](offset + 2),
                **cycle(offset + 3, ";D2;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'PdX2 bisphosphine and acetate {modifier}' : {
            'Reagent': [f'[Pd;D2,D4;+0](-[{fg.X}])-[{fg.X}]', '[O]-[C](=[O])-[C]', fg.bisphosphine] + global_reagent,
            'Exclude_reagent': ['[O;H2]', fg.amine] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                offset + 2 : {
                    'Templates': [
                        f"[O;-1:1]-[C:2](=[O:3])-[C:4].[Pd:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4]",
                         f"[{fg.AM}:7]-[O:1]-[C:2](=[O:3])-[C:4].[Pd:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4].[{fg.AM};+1:7]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Halide exchange',
                },
                offset + 3 : {
                    'Templates': [
                        f"[O;-1:1]-[C:2](=[O:3])-[C:4].[Pd:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4]",
                         f"[{fg.AM}:7]-[O:1]-[C:2](=[O:3])-[C:4].[Pd:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4].[{fg.AM};+1:7]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Halide exchange',
                },
                **catalyst_preparation['Pd(OAc)2 without water']['Bisphosphine'](offset + 4, fg.p_ligand),
                **cycle(offset + 8, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'Pd(OAc)2 monophosphine {modifier}' : {
            'Reagent': [fg.pd_oac_2,f'[{fg.p_ligand}]'] + global_reagent,
            'Exclude_reagent': [fg.bisphosphine,] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 without water']['Monophosphine'](offset + 2, fg.p_ligand),
                offset + 7: {
                    'Templates': [steps["Ligand dissociation"]["Third"]("Pd", fg.p_ligand),],
                    'pKa': [None],
                    'Description': 'Dissociate ligand'
                },
                **cycle(offset + 8, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            }
        },
        f'Pd(OAc)2 bisphosphine {modifier}' : {
            'Reagent': [fg.pd_oac_2, fg.bisphosphine] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                **catalyst_preparation['Attach ligands'](offset, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 without water']['Bisphosphine'](offset + 2, fg.p_ligand),
                **cycle(offset + 6, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            },
        },
        f'Pd-allyl {modifier}' : {
            'Reagent': [f'[{fg.X}]-[Pd;D2,D3]-[#6]-[#6]=[#6]'] + global_reagent,
            'Exclude_reagent': [] + global_exclude,
            'Stages': {
                offset : {
                    'Templates': [steps["Ligand association"]["First"]("Pd", fg.p_ligand, 2),],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                **catalyst_preparation['Pd-allyl'](offset+1),
                **cycle(offset + 3, ";D1;!$([Pd]-[Cl,Br,I,O])"),
            },
        }

    }


catalytic_cycles = {
    "Suzuki" : {
        'Borate' : lambda offset, l = ";D2": 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, None , offset + 10], "Pd;!$([Pd]-[F,Cl,Br,I,O$([O]-S),$([O]-C(=O)C)])", "#5;-1", "[#8:4]", "[#6:5]", "#6", "#6", l1=l), 
                offset + 1: {
                    'Templates' : [steps['Fluoride leaving']['First']],
                    'pKa' : [None],
                    'Description' : 'BF3 to BF2'
                },
                offset + 2: {
                    'Templates': ["[#5;$([#5](-[#9])(-[#9])):2].[OH2:1]>>[#5;$([#5](-[#9])(-[#9])):2].[O;H1;-1:1]"],
                    'pKa': [{'B': 0}],
                    'Description': 'Hydroxide from Water'
                },
                offset + 3 : {
                    'Templates': ["[#5;$([#5](-[#9])(-[#9])):2]-[#9:3].[O;H1;-1:1]>>[#5:2]-[O;+0;H1:1].[#9;-1:3]"],
                    'pKa' : [None],
                    'Description' : 'BF2 to BFOH'
                },
                offset + 4:{
                    'Templates': ["[#5;$([#5](-[#9])(-[#8])):2].[OH2:1]>>[#5:2].[OH1;-1:1]"],
                    'pKa': [{'B': 0}],
                    'Description': 'Hydroxide from Water'
                },
                offset + 5 : {
                    'Templates': ["[#5;$([#5](-[#9])(-[#8])):2]-[#9:3].[O;H1;-1:1]>>[#5:2]-[O;+0;H1:1].[#9;-1:3]"],
                    'pKa' : [None],
                    'Description' : 'BFOH to B(OH)2'
                },
                offset + 6: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]", "[OH1:1]-[#6:2]>>[OH0;-1:1]-[#6:2]"],
                    'pKa': [{'B': -2}, {'B': 0}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide/alkoxide formation' 
                },
                offset + 7: {
                    'Templates': [
                        val for val in steps["Borate formation"].values()
                    ],
                    'pKa': [None, None, None],
                    'Description': 'Borate formation',
                },
                offset + 8: {
                    'Templates': [
                        steps["PdOB formation"]['Borate']
                    ],
                    'pKa': [None],
                    'Description': 'PdOB formation',
                },
                offset + 9: {
                    'Templates' : ["[#6:3]-[Pd:1]-[#8;+1:2]-[#5;-1:4]-[#6;+0:5]>>[#8;+0:2]-[#5;+0:4].[#6:3]-[Pd:1]-[#6:5]"],
                    'pKa': [None],
                    'Description': 'Transmetalation',
                }
            }.items())),

        'Pd-hydroxy' : lambda offset, l = ";D2": 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, None , offset + 10], "Pd;!$([Pd]-[F,Cl,Br,I,O$([O]-S),$([O]-C(=O)C)])", "#5;-1", "[#8:4]", "[#6:5]", "#6", "#6", l1=l),
                offset + 1: {
                    'Templates' : [steps['Fluoride leaving']['First']],
                    'pKa' : [None],
                    'Description' : 'BF3 to BF2'
                },
                offset + 2: {
                    'Templates': ["[#5;$([#5](-[#9])(-[#9])):2].[OH2:1]>>[#5;$([#5](-[#9])(-[#9])):2].[O;H1;-1:1]"],
                    'pKa': [{'B': 0}],
                    'Description': 'Hydroxide from Water'
                },
                offset + 3 : {
                    'Templates': ["[#5;$([#5](-[#9])(-[#9])):2]-[#9:3].[O;H1;-1:1]>>[#5:2]-[O;+0;H1:1].[#9;-1:3]"],
                    'pKa' : [None],
                    'Description' : 'BF2 to BFOH'
                },
                offset + 4:{
                    'Templates': ["[#5;$([#5](-[#9])(-[#8])):2].[OH2:1]>>[#5:2].[OH1;-1:1]"],
                    'pKa': [{'B': 0}],
                    'Description': 'Hydroxide from Water'
                },
                offset + 5 : {
                    'Templates': ["[#5;$([#5](-[#9])(-[#8])):2]-[#9:3].[O;H1;-1:1]>>[#5:2]-[O;+0;H1:1].[#9;-1:3]"],
                    'pKa' : [None],
                    'Description' : 'BFOH to B(OH)2'
                },
                offset + 6: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]", "[OH1:1]-[#6:2]>>[OH0;-1:1]-[#6:2]"],
                    'pKa': [{'B': -2}, {'B': 0}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide/alkoxide formation' 
                },
                offset + 7: {
                    'Templates': [
                        val("Pd") for val in steps["Halide exchange"].values()
                    ],
                    'pKa': [None, None, None, None],
                    'Description': 'Halide exhange',
                },
                offset + 8: {
                    'Templates': [
                        steps["PdOB formation"]['Pd-hydroxy']
                    ],
                    'pKa': [None],
                    'Description': 'PdOB formation',
                },
                offset + 9: {
                    'Templates' : ["[#6:3]-[Pd:1]-[#8;+1:2]-[#5;-1:4]-[#6;+0:5]>>[#8;+0:2]-[#5;+0:4].[#6:3]-[Pd:1]-[#6:5]"],
                    'pKa': [None],
                    'Description': 'Transmetalation',
                }
            }.items())),
    },
    
    "Sonogashira" : lambda offset, l = ";D2":
        dict(sorted({
                **core_cat_cycle([offset, offset + 2, offset + 3], "Pd", "Cu", f"[{fg.X}:4]", "[#6:5]#[#6:6]", "#6", "#6", l1 = l), 
                offset + 1: {
                    'Templates': [
                        f'[Cu:1]-[{fg.X}:2].[#6;H1:3]#[#6;H0:4]>>[Cu:1]-[#6;H0:3]#[#6;H0:4].[{fg.X};H0;-1:2]'
                    ],
                    'pKa': [{'B': -6}],
                    'Description': 'Copper acetylide formation',
                },
                0: {
                'Templates' : [
                    "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
                ],
                'pKa' : [None],
                'Description' : 'CuI2 Decomposition', 
                },
            }.items())),
    "Heck" : lambda offset, l = ";D2":
        {
            offset + 0: {
                'Templates': [ 
                    steps["Oxidative addition"]("Pd", l=l)
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 1: {
                'Templates': [
                    steps["Migratory insertion"]['Normal']("Pd", "#6", "#6", "#6", b_type = "#"),
                    steps["Migratory insertion"]['Normal']("Pd", "#6", "#6", "#6", b_type = "="),
                ],
                'pKa': [None, None],
                'Description': 'Migratory insertion',
            },
            offset + 2: {
                'Templates': [
                    steps["Beta-hydride elimination"]("Pd", "#6", "#6;H2"),
                    steps["Beta-hydride elimination"]("Pd", "#6", "#6;H1")
                ],
                'pKa': [None, None],
                'Description': 'Beta-hydride elimination',
            },
            offset + 3: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride - BA"]("Pd", fg.X)
                ],
                'pKa': [{'B': -6}],
                'Description': 'Reductive elimination',
            },
            offset + 4: {
                'Templates': [ 
                    steps["Oxidative addition"]("Pd", l=l)
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 5: {
                'Templates': [
                    steps["Migratory insertion"]['Normal']("Pd", "#6", "#6", "#6", b_type = "#"),
                    steps["Migratory insertion"]['Normal']("Pd", "#6", "#6", "#6", b_type = "="),
                ],
                'pKa': [None, None],
                'Description': 'Migratory insertion',
            },
            offset + 6: {
                'Templates': [
                    steps["Beta-hydride elimination"]("Pd", "#6", "#6;H2"),
                    steps["Beta-hydride elimination"]("Pd", "#6", "#6;H1")
                ],
                'pKa': [None, None],
                'Description': 'Beta-hydride elimination',
            },
            offset + 7: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride - BA"]("Pd", fg.X)
                ],
                'pKa': [{'B': -6}],
                'Description': 'Reductive elimination',
            },
        },
    "Ullman C-N" : lambda offset, R , l = 3:
        {
            offset: {
                'Templates': [ 
                        f'[#7;D3;H2;+0:3].[Cu;D{l}:1]-[{fg.X}:2]>>[Cu:1]-[#7;H1;+0:3].[{fg.X};H0;-1:2]',
                        f'[#7;D3;H1;+0:3].[Cu;D{l}:1]-[{fg.X}:2]>>[Cu:1]-[#7;H0;+0:3].[{fg.X};H0;-1:2]',
                        f'[#7;D3;H3;+0:3].[Cu;D{l}:1]-[{fg.X}:2]>>[Cu:1]-[#7;H2;+0:3].[{fg.X};H0;-1:2]',
                    ],
                'pKa': [{'B': -6}, {'B': -6}, {'B': -6}],
                'Description': 'Deprotonation',
            },
            offset + 1: {
                'Templates' : [
                    steps["Oxidative addition"]("Cu", f";D{l};!$([Cu]-[{fg.X}])")
                ],
                'pKa' : [None],
                'Description' : 'Oxidative addition', 
            },
            offset + 2: {
                'Templates' : [
                    steps["Reductive elimination"]["Normal"]("Cu", "#6", R)
                ],
                'pKa' : [None],
                'Description' : 'Reductive elimination', 
            }
        },
    "Ullman C-O" : lambda offset, R , l = 3:
        {
            offset : {
                'Templates' : [
                    '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                    #'[n;H1:1]:[n;H0:2]>>[n;H1:1]:[n;H0:2]'
                ],
                'pKa' : [None],
                'Description' : 'Tautomerization'
            },
            offset + 1: {
                'Templates' : [
                    '[O;H1:1]>>[O;H0;-1:1]'
                ],
                'pKa' : [{'B' : 3}],
                'Description' : 'Deprotonation'
            },
            offset + 2: {
                'Templates' : [
                    f'[O;-1:1].[Cu;D{l}:2]-[{fg.X}:3]>>[O;+0:1]-[Cu:2].[{fg.X};-1:3]',
                    f'[O:1]-[{fg.AM}:4].[Cu;D{l}:2]-[{fg.X}:3]>>[O;+0:1]-[Cu:2].[{fg.X};-1:3].[{fg.AM};+1:4]'
                ],
                'pKa' : [None, None],
                'Description' : 'Coordination of alkoxide'
            },
            offset + 3: {
                'Templates' : [
                    steps["Oxidative addition"]("Cu", f";D{l};!$([Cu]-[{fg.X}])")
                ],
                'pKa' : [None],
                'Description' : 'Oxidative addition', 
            },
            offset  + 4: {
                'Templates' : [
                    steps["Reductive elimination"]["Normal"]("Cu", "#6", R)
                ],
                'pKa' : [None],
                'Description' : 'Reductive elimination', 
            }
        },
    "Migita thioether" : lambda offset, l = ";D2" : {
        **core_cat_cycle([offset, None, offset + 3], "Pd", "Cu", None, None, "#6", "#16", l1 = l), 
        offset + 1 : {
                'Templates' : [
                    '[S;H1:1]>>[S;H0;-1:1]'
                ],
                'pKa' : [{'B' : 0}],
                'Description' : 'Deprotonation'
        },
        offset + 2: {
            'Templates' : [
                f'[S;-1:1].[Pd:2]-[{fg.X}:3]>>[S;+0:1]-[Pd:2].[{fg.X};-1:3]',
                f'[S:1]-[{fg.AM}:4].[Pd:2]-[{fg.X}:3]>>[S;+0:1]-[Pd:2].[{fg.X};-1:3].[{fg.AM};+1:4]'
            ],
            'pKa' : [None, None],
            'Description' : 'Coordination of thiolate'
        },

    },
    "Alkoxycarbonylation" : lambda offset, l = ";D2" : {
        **core_cat_cycle([offset, None, offset + 4], "Pd", "Cu", None, None, "#6", "#6", l1 = l), 
        offset + 1 : {
                'Templates' : [
                    '[O;H1:1]>>[O;H0;-1:1]'
                ],
                'pKa' : [{'B' : 0}],
                'Description' : 'Deprotonation'
        },
        offset + 2: {
            'Templates' : [
                f'[O;-1:1].[Pd:2]-[{fg.X}:3]>>[O;+0:1]-[Pd:2].[{fg.X};-1:3]',
                f'[O:1]-[{fg.AM}:4].[Pd:2]-[{fg.X}:3]>>[O;+0:1]-[Pd:2].[{fg.X};-1:3].[{fg.AM};+1:4]'
            ],
            'pKa' : [None, None],
            'Description' : 'Coordination of alkoxide'
        },
        offset + 3: {
            'Templates' : [
                f'[Pd:1]-[O:2].[C:3]#[O:4]>>[Pd:1]-[C;+0:3](=[O;+0:4])-[O:2]'
            ],
            'pKa' : [None],
            'Description' : '1,1 migratory insertion'
        }
    },
    "Aminocarbonylation" : lambda offset, l = ";D2" : {
        **core_cat_cycle([offset, None, offset + 3], "Pd", "Cu", None, None, "#6", "#6", l1 = l), 
        offset + 1 : {
            'Templates': [
                f'[#7;H1:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H0:1]-[Pd:2]',
                f'[#7;H2:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H1:1]-[Pd:2]',
                f'[#7;H3:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H2:1]-[Pd:2]',
            ],
            'pKa': [{'B': 0}, {'B': 0}, {'B': 0}], #??
            'Description': 'Deprotonation',
        },
        offset + 2: {
            'Templates' : [
                f'[Pd:1]-[#7:2].[C:3]#[O:4]>>[Pd:1]-[C;+0:3](=[O;+0:4])-[#7:2]'
            ],
            'pKa' : [None],
            'Description' : '1,1 migratory insertion'
        }
    },

    "Alkyne hydrogenation" : lambda offset, M : {
            offset : {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 1: {
                'Templates': [
                    steps["Migratory insertion"]['Hydride'](M, 2, f"#6;H{num_H}", "#6", b_type = "#") for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Migratory insertion',
            },
            offset + 2: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride"](M, f"#6;H{num_H}")  for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Reductive elimination',
            },
            offset + 3: {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 4: {
                'Templates': [
                    steps["Migratory insertion"]['Hydride'](M, 2, f"#6;H{num_H}", "#6", b_type = "=") for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Migratory insertion',
            },
            offset + 5: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride"](M, f"#6;H{num_H}")  for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Reductive elimination',
            },
        },
    "Alkene hydrogenation" : lambda offset, M : {
            offset : {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 1: {
                'Templates': [
                    steps["Migratory insertion"]['Hydride'](M, 2, f"#6;H{num_H}", "#6", b_type = "=") for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Migratory insertion',
            },
            offset + 2: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride"](M, f"#6;H{num_H}")  for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Reductive elimination',
            },
        },
    "Debenzylation" : lambda offset : {
            offset : {
                'Templates': [f'[Pd;D0;+0:1].[H:2]-[H:3]>>[H:2]-[Pd:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
            },
            offset + 1: {
                'Templates': [
                    '[O,N:1]-[C:2]-[c:3].[H:4]-[Pd:5]-[H:6]>>[H:4]-[Pd:5]-[C:2]-[c:3].[O,N:1]-[H:6]'
                ],
                'pKa': [None],
                'Description': 'Hydrogenolysis',
            },
            offset + 2: {
                'Templates': [
                    steps["Reductive elimination"]["Hydride"]("Pd", f"#6;H{num_H}")  for num_H in range(3)
                ],
                'pKa': [None, None, None],
                'Description': 'Reductive elimination',
            },
        },
    "Buchwald-Hartwig amination" : lambda offset, l = ";D1" : 
        dict(sorted({ # default one ligand on Pd
            **core_cat_cycle([offset, None , offset + 2], "Pd;!$([Pd]-[F,Cl,Br,I,O$([O]-S),$([O]-C(=O)C)])", "", "", "", "#7", "#6", l1=l),
            offset + 1 : {
                'Templates': [
                    f'[#7;H1:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H0:1]-[Pd:2]',
                    f'[#7;H2:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H1:1]-[Pd:2]',
                    f'[#7;H3:1].[Pd:2]-[{fg.X}:3]>>[{fg.X};-1:3].[#7;H2:1]-[Pd:2]',
                ],
                'pKa': [{'B': 0}, {'B': 0}, {'B': 0}], #??
                'Description': 'Deprotonation',
            },
        }.items())),
    "Miyaura borylation" : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, offset + 2 , offset + 3], "Pd;!$([Pd]-[F,Cl,Br,I,O$([O]-S)])", "#5", "[C:4]-[C:5](=[O:6])-[O:7]", "[#5:8]", "#5", "#6", l1=l),
                offset + 1 : {
                    'Templates': [
                        f"[O;-1:1]-[C:2](=[O:3])-[C:4].[Pd{more_subs(l)}:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4]",
                         f"[{fg.AM}:7]-[O:1]-[C:2](=[O:3])-[C:4].[Pd{more_subs(l)}:5]-[{fg.X}:6]>>[{fg.X};-1:6].[Pd:5]-[O;+0:1]-[C:2](=[O:3])-[C:4].[{fg.AM};+1:7]"
                    ],
                    'pKa': [None, None], #??
                    'Description': 'Halide exchange',
                }
        }.items())),
    "Nitro to amine" : lambda offset, M : {
        offset : {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
        },
        offset + 1: {
            'Templates': [ f"[N;+1:1](=[O:2])-[O;-1:3].[H:4]-[{M}:5]-[H:6]>>[H:4]-[N;+1:1](-[O:2]-[{M}:5]-[H:6])-[O;-1:3]"
            ],
            'pKa': [None],
            'Description': 'Migratory insertion',
        },
        offset + 2: {
            'Templates': [
                steps["Reductive elimination"]["Hydride"](M, f"#8;H0")
            ],
            'pKa': [None],
            'Description': 'Reductive elimination',
        },
        offset + 3:{
            'Templates' : [
                "[N;H1;+1:1](-[O;H1:2])-[O;-1:3]>>[N;H0;+0:1](-[O;H1:2])-[O;+0;H1:3]"
            ],
            'pKa': [None],
            'Description': 'Proton transfer',
        },
        offset + 4:{
            'Templates' : [
                "[N;H0;+0:1](-[O;H1:2])-[O;+0;H1:3]>>[N;H0;+0:1](-[O;H0;-1:2])-[O;+1;H2:3]"
            ],
            'pKa': [None],
            'Description': 'Proton transfer',
        },
        offset + 5:{
            'Templates' : [
                "[N;H0;+0:1](-[O;H0;-1:2])-[O;+1;H2:3]>>[O;H2;+0:3].[N:1]=[O;+0:2]"
            ],
            'pKa': [None],
            'Description': 'N=O forms',
        },
        offset + 6: {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
        },
        offset + 7: {
            'Templates': [ f"[N:1]=[O:2].[H:4]-[{M}:5]-[H:6]>>[H:4]-[N:1]-[O:2]-[{M}:5]-[H:6]"
            ],
            'pKa': [None],
            'Description': 'Migratory insertion',
        },
        offset + 8: {
            'Templates': [
                steps["Reductive elimination"]["Hydride"](M, f"#8;H0")
            ],
            'pKa': [None],
            'Description': 'Reductive elimination',
        },
        offset + 9: {
                'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
                ],
                'pKa': [None],
                'Description': 'Oxidative addition',
        },
        offset + 10:{
            'Templates' : [
                f"[N;H1:1]-[O;H1:2].[H:4]-[{M}:5]-[H:6]>>[N;H1:1]-[H:4].[O;H1:2]-[{M}:5]-[H:6]"
            ],
            'pKa': [None],
            'Description': 'Hydrogenolysis',
        },
        offset + 11: {
            'Templates': [
                steps["Reductive elimination"]["Hydride"](M, f"#8;H1")
            ],
            'pKa': [None],
            'Description': 'Reductive elimination',
        },
    },
    "Azide to amine" : lambda offset, M : {
        offset : {
            'Templates': [f'[{M};D0;+0:1].[H:2]-[H:3]>>[H:2]-[{M}:1]-[H:3]'
            ],
            'pKa': [None],
            'Description': 'Oxidative addition',
        },
        offset + 1 : {
            'Templates': [f'[H:2]-[{M}:1]-[H:3].[N:4]=[N;+1:5]=[N;-1:6]>>[H:2]-[N:4]-[N;+1:5](=[N;-1:6])-[{M}:1]-[H:3]'
            ],
            'pKa': [None],
            'Description': 'Migratory insertion',
        },
        offset + 2 : {
            'Templates': [f"[N;+1:5]-[{M}:1]-[H:3]>>[{M}:1].[N;+1:5]-[H:3]"],
            'pKa': [None],
            'Description': 'Reductive elimination',
        },
        offset + 3 : {
            'Templates': ['[H:1]-[N:2]-[N;+1:3]-[H:4]>>[H:1]-[N;+1:2](-[H:4])-[N;+0:3]'
            ],
            'pKa': [None],
            'Description': 'Proton shift',
        },
        offset + 4 : {
            'Templates': ['[N;+1:2]-[N:3]=[N;-1:4]>>[N:3]#[N;+0:4].[N;+0:2]'
            ],
            'pKa': [None],
            'Description': 'N2 leaves',
        }
    },
    'CMD C-H arylation' : lambda offset, l: {
        offset : {
            'Templates' : [steps["Oxidative addition"]("Pd", l=l)],
            'pKa' : [None],
            'Description' : 'Oxidative addition',
        },
        offset + 1 : {
            'Templates' : [f"[Pd:1]-[{fg.X}:2].[O;-1:3]-[C:4]=[O:5]>>[Pd:1]-[O;+0:3]-[C:4]=[O:5].[{fg.X};-1:2]",
                           f"[Pd:1]-[{fg.X}:2].[{fg.AM}:6]-[O:3]-[C:4]=[O:5]>>[{fg.X};-1:2].[{fg.AM};+1:6].[Pd:1]-[O:3]-[C:4]=[O:5]",
                           f"[Pd:1]-[{fg.X}:2].[Ag:6]-[O:3]-[C:4]=[O:5]>>[Ag:6]-[{fg.X}:2].[Pd:1]-[O:3]-[C:4]=[O:5]"],
            'pKa' : [None, None, None],
            'Description' : 'Halide exchange',
        },
        offset + 2 : {
            'Templates' : ["[Pd:1]-[O;+0:2]-[C:3]=[O:4].[c:5]-[H:6]>>[O:2]=[C:3]-[O:4]-[H:6].[c:5]-[Pd:1]"],
            'pKa' : [None],
            'Description' : "Concerted metalation deprotonation"
        },
        offset + 3 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Pd", "#6", "#6")],
            'pKa' : [None],
            'Description' : "Reductive elimination"
        }
    },
    "Chan-Lam C-N" : lambda offset : {
        offset : {
            'Templates' : ["[Cu:1]-[O;+0:2]-[C:3]=[O:4].[#7:5]-[H:6]>>[O:2]=[C:3]-[O:4]-[H:6].[#7:5]-[Cu:1]"],
            'pKa' : [None],
            'Description' : "Exchange with amine"
        },
        offset + 1: {
            'Templates' : [steps['Transmetalation']("Cu", "#5", "[#8:4]", "[#6:5]", S="#7")],
            'pKa' : [None], 
            'Description' : "Transmetalation"
        },
        offset + 2 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Cu", "#7", "#6")],
            'pKa' : [None],
            'Description' : "Reductive elimination"
        }
    },
    "Chan-Lam C-O" : lambda offset : {
        offset : {
            'Templates' : ["[Cu:1]-[O;+0:2]-[C:3]=[O:4].[#8:5]-[H:6]>>[O:2]=[C:3]-[O:4]-[H:6].[#8:5]-[Cu:1]"],
            'pKa' : [None],
            'Description' : "Exchange with amine"
        },
        offset + 1: {
            'Templates' : [steps['Transmetalation']("Cu", "#5", "[#8;$([#8]-[#6]=[#8]):4]", "[#6:5]", S="#8")],
            'pKa' : [None], 
            'Description' : "Transmetalation"
        },
        offset + 2 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Cu", "#8", "#6")],
            'pKa' : [None],
            'Description' : "Reductive elimination"
        }
    },
    "Rosenmund" : lambda offset : {
        offset : {
            'Templates' : [f"[Cu:1]-[C:2]#[N:3].[#6:4]-[{fg.X}:5]>>[Cu:1](-[C:2]#[N:3])(-[{fg.X}:5])-[#6:4]"],
            'pKa' : [None], 
            'Description' : "Oxidative addition"
        },
        offset + 1 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Cu", "#6;$([#6]#[N])", "#6")],
            'pKa' : [None],
            'Description' : ''
        }
    },
    'Hiyama' : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, offset + 2 , offset + 3], "Pd", "#14;-1", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
                offset + 1 : {
                    'Templates': [
                        "[Si:1].[F;-1:2]>>[Si;-1:1]-[F;+0:2]",
                        f"[Si:1].[F:2]-[{fg.AM}:3]>>[Si;-1:1]-[F;+0:2].[{fg.AM};+1:3]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Silane activation',
                }
        }.items())),
    'Hiyama with copper' : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, offset + 3 , offset + 4], "Pd", "Cu", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
                offset + 1 : {
                    'Templates': [
                        "[Si:1].[F;-1:2]>>[Si;-1:1]-[F;+0:2]",
                        f"[Si:1].[F:2]-[{fg.AM}:3]>>[Si;-1:1]-[F;+0:2].[{fg.AM};+1:3]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Silane activation',
                },
                offset + 2 : {
                    'Templates' : [f"[Cu:1]-[{fg.X}:2].[Si;-1:3]-[#6:4]>>[Cu:1]-[#6:4].[Si;-1:3]-[{fg.X}:2]"],
                    'pKa': [None],
                    'Description': 'Si to Cu transmetalation',
                }
        }.items())),
    'Hiyama with copper, electron rich Si' : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, offset + 2, offset + 3], "Pd", "Cu", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
                offset + 1 : {
                    'Templates' : [f"[Cu:1]-[{fg.X}:2].[Si:3](-[O:7])(-[O:8])(-[O:9])-[#6:4]>>[Cu:1]-[#6:4].[Si:3](-[O:7])(-[O:8])(-[O:9])-[{fg.X}:2]"],
                    'pKa': [None],
                    'Description': 'Si to Cu transmetalation',
                }
        }.items())),
    'Hiyama base activation' : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, None , offset + 5], "Pd", "#14;-1", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
                offset + 1: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]", "[OH1:1]-[#6:2]>>[OH0;-1:1]-[#6:2]"],
                    'pKa': [{'B': -2}, {'B': 0}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide/alkoxide formation' 
                },
                offset + 2 : {
                    'Templates': [
                        "[Si:1].[O;-1:2]>>[Si;-1:1]-[O;+0:2]",
                        f"[Si:1].[O:2]-[{fg.AM}:3]>>[Si;-1:1]-[O;+0:2].[{fg.AM};+1:3]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Silane activation',
                },
                offset + 3 : {
                    'Templates' : [f"[Si;-1:1]-[O:2].[Pd:3]-[{fg.X}:4]>>[Si;-1:1]-[O;+1:2]-[Pd:3].[{fg.X};-1:4]"],
                    'pKa': [None],
                    'Description': 'Si-O-Pd formation',
                },
                offset + 4 : {
                    'Templates' : [f"[#6:4]-[Si;-1:1]-[O;+1:2]-[Pd:3]>>[Pd:3]-[#6:4].[Si;+0:1]-[O;+0:2]"],
                    'pKa': [None],
                    'Description': 'Transmetalation',
                }
        }.items())),
    'Hiyama base activation with Cu' : lambda offset, l = ";D2" : 
        dict(sorted({ # default two ligands on Pd
                **core_cat_cycle([offset, offset + 5 , offset + 6], "Pd", "Cu", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
                offset + 1: {
                    'Templates': ["[OH2:1]>>[OH1;-1:1]", "[OH1:1]-[#6:2]>>[OH0;-1:1]-[#6:2]"],
                    'pKa': [{'B': -2}, {'B': 0}], # might have to adjust this (this is kinda arbitrary)
                    'Description': 'Hydroxide/alkoxide formation' 
                },
                offset + 2 : {
                    'Templates': [
                        "[Si:1].[O;-1:2]>>[Si;-1:1]-[O;+0:2]",
                        f"[Si:1].[O:2]-[{fg.AM}:3]>>[Si;-1:1]-[O;+0:2].[{fg.AM};+1:3]"
                    ],
                    'pKa': [None, None],
                    'Description': 'Silane activation',
                },
                offset + 3 : {
                    'Templates' : [f"[Si;-1:1]-[O:2].[Cu:3]-[{fg.X}:4]>>[Si;-1:1]-[O;+1:2]-[Cu:3].[{fg.X};-1:4]"],
                    'pKa': [None],
                    'Description': 'Si-O-Cu formation',
                },
                offset + 4 : {
                    'Templates' : [f"[#6:4]-[Si;-1:1]-[O;+1:2]-[Cu:3]>>[Cu:3]-[#6:4].[Si;+0:1]-[O;+0:2]"],
                    'pKa': [None],
                    'Description': 'Transmetalation',
                }
        }.items())),
    'Pd alcohol dehalogenation' : lambda offset, l = ";D2" : {
        offset : {
            'Templates' : [steps["Oxidative addition"]("Pd", l = l)], 
            'pKa' : [None],
            'Description' : 'Oxidative addition'
        }, 
        offset + 1 : {
            'Templates' : ["[O;H1:1]>>[O;H0;-1:1]"],
            'pKa' : [{'B' : -2}],
            'Description' : 'Deprotonation'
        },
        offset + 2 : {
            'Templates' : [steps["Halide exchange"]["Alkoxide"]("Pd"),
                            steps["Halide exchange"]["Alkoxide-Metal"]("Pd"),],
            'pKa' : [None, None],
            'Description' : 'Halide exchange'
        },
        offset + 3 : {
            'Templates' : [steps["Beta-hydride elimination"]("Pd", "#8", f"C;H{n+1}") for n in range(3)],
            'pKa' : [None, None, None],
            'Description' : 'Beta-hydride elimination',
        },
        offset + 4 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Pd", "#6", "H")],
            'pKa' : [None],
            'Description' : 'Reductive elimination'
        }   
    },
    'Pd dehalogenation with hydrogen' : lambda offset, l = ';D2' : {
        offset : {
            'Templates' : [steps["Oxidative addition"]("Pd", l = l)], 
            'pKa' : [None],
            'Description' : 'Oxidative addition'
        }, 
        offset + 1 : {
            'Templates': [f'[Pd;D0;+0:1].[H:2]-[H:3]>>[H:2]-[Pd:1]-[H:3]',
                          f'[Pd;D0;+0:1].[Si:2]-[H:3]>>[Si:2]-[Pd:1]-[H:3]'
                ],
            'pKa': [None, None],
            'Description': 'Oxidative addition',
        },
        offset + 2 : {
            'Templates' : [steps['Transmetalation']("Pd;D2", "Pd;D2", f"[{fg.X}:6]", "[H:5]")],
            'pKa' : [None],
            'Description' : 'Transmetalation'
        },
        offset + 3 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Pd", "#6", "H")],
            'pKa' : [None],
            'Description' : 'Reductive elimination'
        },
        offset + 4: {
            'Templates': [steps["Reductive elimination"]["Hydride - BA"]("Pd", fg.X)],
            'pKa': [{'B': -6}],
            'Description': 'Reductive elimination',
        },
    },
    'Pd formate hydrogenation' : lambda offset, l = ";D2" : {
        offset : {
            'Templates' : [steps["Oxidative addition"]("Pd", l = l)], 
            'pKa' : [None],
            'Description' : 'Oxidative addition'
        }, 
        offset + 1 : {
            'Templates' : ["[C:2](=[O:3])-[O;H1:1]>>[C:2](=[O:3])-[O;H0;-1:1]"],
            'pKa' : [{'B' : -2}],
            'Description' : 'Deprotonation'
        },
        offset + 2 : {
            'Templates' : [f"[Pd:1]-[{fg.X}:2].[C:3](=[O:4])-[O;H0;-1:5]>>[C:3](=[O:4])-[O;H0;+0:5]-[Pd:1].[{fg.X};-1:2]",
                           f"[Pd:1]-[{fg.X}:2].[C:3](=[O:4])-[O;H0:5]-[{fg.AM}:6]>>[C:3](=[O:4])-[O;H0;+0:5]-[Pd:1].[{fg.X};-1:2].[{fg.AM};+1:6]"],
            'pKa' : [None, None],
            'Description' : 'Halide exchange'
        },
        offset + 3: {
            'Templates' : ["[Pd:1]-[O:2]-[C:3](=[O:4])-[H:5]>>[Pd:1]-[H:5].[O:2]=[C:3]=[O:4]"],
            'pKa' : [None],
            'Description' : 'Decarboxylation'
        },
        offset + 4 : {
            'Templates' : [steps["Reductive elimination"]["Normal"]("Pd", "#6", "H")],
            'pKa' : [None],
            'Description' : 'Reductive elimination'
        }   
    },

    'Alpha arylation' : lambda offset, l = ";D2" : {
        **core_cat_cycle([offset, None , offset + 4], "Pd", "#14;-1", f"[{fg.X}:5]", "[#6:6]", "#6", "#6", l1=l),
        offset + 1 : {
            'Templates': ["[O:1]=[C:2]-[C;H1:3]>>[O;-1:1]-[C:2]=[C;H0:3]",
                        "[O:1]=[C:2]-[C;H2:3]>>[O;-1:1]-[C:2]=[C;H1:3]",
                        "[O:1]=[C:2]-[C;H3:3]>>[O;-1:1]-[C:2]=[C;H2:3]"
                            ],
            'pKa': [{'B' : 3},{'B' : 3},{'B' : 3}],
            'Description': 'Deprotonate'
        },
        offset + 2 : {
            'Templates' : [f"[Pd:1]-[{fg.X}:2].[C:3](=[C:4])-[O;H0;-1:5]>>[C:3](=[C:4])-[O;H0;+0:5]-[Pd:1].[{fg.X};-1:2]"],
            'pKa' : [None],
            'Description' : 'Halide exchange'
        },
        offset + 3 : {
            'Templates' : ["[C:3](=[C:4])-[O;H0;+0:5]-[Pd:1]>>[Pd:1]-[C:4]-[C:3]=[O:5]"],
            'pKa' : [None],
            'Description' : 'Beta-carbon elimination'
        }
    }

}

class_reaction_templates = {
    (
        'Bromo Suzuki coupling',
        'Chloro Suzuki coupling',
        'Iodo Suzuki coupling',
        'Triflyloxy Suzuki coupling',
        'Tosyloxy Suzuki coupling',
        'Mesyloxy Suzuki coupling',
        'Triflyl Suzuki coupling',
    ): 
    {
        'Pd(0) ligandless - borate': {
            'Reagent': ['[Pd;D0;+0]',fg.oh_or,],
            'Exclude_reagent': [f'[{fg.p_ligand}]'],
            'Stages': catalytic_cycles['Suzuki']['Borate'](0, ";D0"),
        },
        'Pd(0) ligandless - hydroxy': {
            'Reagent': ['[Pd;D0;+0]',fg.oh_or,],
            'Exclude_reagent': [f'[{fg.p_ligand}]'],
            'Stages': catalytic_cycles['Suzuki']['Pd-hydroxy'](0, ";D0"),
        },
        'Pd(0) with ligand - borate': {
            'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])"),
                **catalytic_cycles['Suzuki']['Borate'](2, ";D2")},
        },
        'Pd(0) with ligand - hydroxy': {
            'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])"),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](2, ";D2")},
        },
        'Pd(0) with single ligand - borate': {
            'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                0 : {
                    'Templates': [steps["Ligand association"]["First"]("Pd", fg.p_ligand, 0, "!$([P]-[Pd])"),],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                **catalytic_cycles['Suzuki']['Borate'](1, ";D1")},
        },
        'Pd(0) with single ligand - hydroxy': {
            'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                0 : {
                    'Templates': [steps["Ligand association"]["First"]("Pd", fg.p_ligand, 0, "!$([P]-[Pd])"),],
                    'pKa': [None],
                    'Description': 'Attach first ligand'
                },
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](1, ";D1")},
        },
        'Pd(0) with three or four ligands - borate': {
            'Reagent': [f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](2, ";D2")},
        },
        'Pd(0) with three or four ligands - hydroxy': {
            'Reagent': [f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](2, ";D2")},
        },
        'Pd(0) with two ligands - borate': {
            'Reagent': [f'[Pd;D2;$([Pd](-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalytic_cycles['Suzuki']['Borate'](0, ";!$([Pd]-[F,Cl,Br,I,O$([O]-S),$([O]-C(=O)C)]);D2")},
        },
        'Pd(0) with two ligands - hydroxy': {
            'Reagent': [f'[Pd;D2;$([Pd](-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](0, ";!$([Pd]-[F,Cl,Br,I,O$([O]-S),$([O]-C(=O)C)]);D2")},
        },
        'Pd(0) with two, three, or four ligands - borate' : {
            'Reagent': [f'[Pd;D2,D3,D4;$([Pd](-[*])-[*]),$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Remove three ligands'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](3, ";D1")},
        },
        'Pd(0) with two, three, or four ligands - hydroxy' : {
            'Reagent': [f'[Pd;D2,D3,D4;$([Pd](-[*])-[*]),$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Remove three ligands'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](3, ";D1")},
        },
        'PdX2 with phosphine, base, water - borate' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', f'[{fg.p_ligand}]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine'](2, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](7, ";D2")},
        },
        'PdX2 with phosphine, base, water - hydroxy' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', f'[{fg.p_ligand}]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                 **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine'](2, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](7, ";D2")},
        },
        'PdX2 with phosphine, 1 ligand, base, water - borate' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', f'[{fg.p_ligand}]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine'](2, "Pd", fg.p_ligand),
                7 : {
                    'Templates': [steps["Ligand dissociation"]["Third"]("Pd", fg.p_ligand),],
                    'pKa': [None],
                    'Description': 'Remove second ligand'
                },
                **catalytic_cycles['Suzuki']['Borate'](8, ";D1")},
        },
        'PdX2 with phosphine, 1 ligand, base, water - hydroxy' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', f'[{fg.p_ligand}]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                 **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine'](2, "Pd", fg.p_ligand),
                7 : {
                    'Templates': [steps["Ligand dissociation"]["Third"]("Pd", fg.p_ligand),],
                    'pKa': [None],
                    'Description': 'Remove second ligand'
                },
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](8, ";D1")},
        },
        'PdX2 with 2 eq phosphine, base, water - borate' : {
            'Reagent': ['[Pd;+0;D4;$([Pd](-[P])(-[P])(-[Cl,Br,I])(-[Cl,Br,I]))]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine - 2 eq'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](5, ";D1")},
        },
        'PdX2 with 2 eq phosphine, base, water - hydroxy' : {
            'Reagent': ['[Pd;+0;D4;$([Pd](-[P])(-[P])(-[Cl,Br,I])(-[Cl,Br,I]))]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['PdX2 with phosphines water and base']['Monophosphine - 2 eq'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](5, ";D1")},
        },
        'PdX2 with bisphosphine, base, water - borate' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', fg.bisphosphine, fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                 **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Bisphosphine'](2, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](9, ";D2")},
        },
        'PdX2 with bisphosphine, base, water - hydroxy' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', fg.bisphosphine, fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                 **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with phosphines water and base']['Bisphosphine'](2, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](9, ";D2")},
        },
        'PdX2 with 1 eq. bisphosphine, base, water - borate' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D4;$([Pd]1-[P]-[*]-[*]-[*]-[P]-1)][Cl,Br,I]', fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['PdX2 with phosphines water and base']['Bisphosphine - 1 eq'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](5, ";D2")},
        },
        'PdX2 with 1 eq. bisphosphine, base, water - hydroxy' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D4;$([Pd]1-[P]-[*]-[*]-[*]-[P]-1)][Cl,Br,I]', fg.oh_or],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['PdX2 with phosphines water and base']['Bisphosphine - 1 eq'](0, "Pd", fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](5, ";D2")},
        },
        'Pd Buchwald precatalyst - borate' : {
            'Reagent': [f'[{fg.anion}]-[Pd](-[P;H0])-[#6]:[#6]-[#6]:[#6]-[#7;H2,H1]', fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Buchwald precatalyst'](0),
                **catalytic_cycles['Suzuki']['Borate'](2, ";D1")},
        },
        'Pd Buchwald precatalyst - hydroxy' : {
            'Reagent': [f'[{fg.anion}]-[Pd](-[P;H0])-[#6]:[#6]-[#6]:[#6]-[#7;H2,H1]',fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Buchwald precatalyst'](0),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](2, ";D1")},
        },
        'Pd amine activation - borate' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', '[C;H3]-[C;H2]-[N](-[C;H2]-[C;H3])-[C;H2]-[C;H3]', f'[{fg.p_ligand}]' ,fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with amine base'](2),
                **catalytic_cycles['Suzuki']['Borate'](5, ";D2")
            },
        },
        'Pd amine activation - hydroxy' : {
            'Reagent': ['[Cl,Br,I][Pd;+0;D2][Cl,Br,I]', '[C;H3]-[C;H2]-[N](-[C;H2]-[C;H3])-[C;H2]-[C;H3]', f'[{fg.p_ligand}]', fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['PdX2 with amine base'](2),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](5, ";D2")
            },
        },
        'Pd(OAc)2 with phosphine - hydroxy' : {
            'Reagent': [fg.pd_oac_2, f'[{fg.p_ligand}]', fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 with phosphines']['Monophosphine'](2, fg.p_ligand),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](8, ";D2")},
        },
        'Pd(OAc)2 with phosphine - borate' : {
            'Reagent': [fg.pd_oac_2, f'[{fg.p_ligand}]' ,fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 with phosphines']['Monophosphine'](2, fg.p_ligand),
                **catalytic_cycles['Suzuki']['Borate'](8, ";D2")},
        },
        'Pd(OAc)2 with bisphosphine - hydroxy' : {
            'Reagent': [fg.pd_oac_2, fg.bisphosphine, fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 with phosphines']['Bisphosphine'](2,),
                **catalytic_cycles['Suzuki']['Pd-hydroxy'](8, ";D2")},
        },
        'Pd(OAc)2 with bisphosphine - borate' : {
            'Reagent': [fg.pd_oac_2,  fg.bisphosphine ,fg.oh_or,],
            'Exclude_reagent': ['[Fe;D0]'],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 with phosphines']['Bisphosphine'](2,),
                **catalytic_cycles['Suzuki']['Borate'](8, ";D2")},
        },
    },

    (
        'Bromo Sonogashira coupling',
        'Iodo Sonogashira coupling',
        'Chloro Sonogashira coupling',
        'Triflyloxy Sonogashira coupling',
    ): 
    {   **general_pd_coupling(catalytic_cycles['Sonogashira'], offset=1, global_reagent=[f'[Cu]-[{fg.X}]'], global_exclude=[fg.nhc,'[Fe;D0]'])

    },
    ( 'Negishi coupling'
    ): {
        **general_pd_coupling(lambda o, l: core_cat_cycle([o, o+1, o+2], "Pd", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = l), offset = 0, global_reagent=['[Zn]~[#6]'], global_exclude=[fg.nhc, '[Fe;D0]']),
        **general_pd_coupling((lambda o, l: 
            {
                0 : {
                    'Templates' : [
                        '[#6:1]-[Br,I:2].[Zn:3]>>[#6:1]-[Zn:3]-[Br,I:2]'
                    ],
                    'pKa' : [None],
                    'Description' : 'In-situ organozinc preparation', 
                },

                # 0 : {
                #     'Templates' : [
                #         steps["Oxidative addition"]("Zn")
                #     ],
                #     'pKa' : [None],
                #     'Description' : 'In-situ organozinc preparation', 
                # },
                **core_cat_cycle([o, o+1, o+2], "Pd", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = l), 
            }),
            offset = 1 , global_reagent=['[Zn;D0]'], global_exclude=[fg.nhc, '[Fe;D0]', '[Zn]~[#6]'], modifier='in-situ')
                #**core_cat_cycle([0, 1, 2], "Pd", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D0")
    },
    ('Palladium Kumada coupling') : {
        'Pd(0) with three or four ligands - borate': {
            'Reagent': [f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]','[Mg]',],
            'Exclude_reagent': ['[Zn]'],
            'Stages': {
                **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
                **core_cat_cycle([2, 3, 4], "Pd", "Mg", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2;!$([Pd]-[Cl,Br,I,O])")
            },
        },
    },
    
    ('Nickel Negishi coupling') : 
    {
        'NiX2P2 and Zn' : {
            'Reagent': [f'[Ni](-[{fg.X}])(-[{fg.X}])',  f'[{fg.p_ligand}]' , '[Zn;D0]',],
            'Exclude_reagent': [fg.amidine_ligand,],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Ni", fg.p_ligand , c = 2),
                2 : {
                    'Templates' : [f'[Zn;D0:1].[Ni;D4:2](-[{fg.X}:3])-[{fg.X}:4]>>[Zn:1](-[{fg.X}:3])-[{fg.X}:4].[Ni:2]'],
                    'pKa' : [None],
                    'Description' : 'Reduction of Ni(II)', 
                },
                3 : {
                    'Templates' : [
                        '[#6:1]-[Br,I:2].[Zn;D0:3]>>[#6:1]-[Zn:3]-[Br,I:2]'
                    ],
                    'pKa' : [None],
                    'Description' : 'In-situ organozinc preparation', 
                },
                **core_cat_cycle([4, 5, 6], "Ni", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2;!$([Ni]-[Cl,Br,I,O])")
            }
        },
        'NiX2N2 and Zn' : {
            'Reagent': [f'[Ni](-[{fg.X}])(-[{fg.X}])',  f'[{fg.n_n_ligand}]' , '[Zn;D0]',],
            'Exclude_reagent': [fg.amidine_ligand,],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Ni", f"#7;{fg.n_n_ligand}", c = 2),
                2 : {
                    'Templates' : [f'[Zn;D0:1].[Ni;D4:2](-[{fg.X}:3])-[{fg.X}:4]>>[Zn:1](-[{fg.X}:3])-[{fg.X}:4].[Ni:2]'],
                    'pKa' : [None],
                    'Description' : 'Reduction of Ni(II)', 
                },
                3 : {
                    'Templates' : [
                        '[#6:1]-[Br,I:2].[Zn;D0:3]>>[#6:1]-[Zn:3]-[Br,I:2]'
                    ],
                    'pKa' : [None],
                    'Description' : 'In-situ organozinc preparation', 
                },
                **core_cat_cycle([4, 5, 6], "Ni", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2;!$([Ni]-[Cl,Br,I,O])")
            }
        },
        'NiX2P2 and organozinc' : {
            'Reagent': [f'[Ni](-[{fg.X}])(-[{fg.X}])',  f'[{fg.p_ligand}]' , '[Zn;D2]',],
            'Exclude_reagent': [fg.amidine_ligand,],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Ni", fg.p_ligand , c = 0),
                2 : {
                    'Templates' : [steps['Transmetalation']( "Ni;D4", "Zn", f"[{fg.X}:4]", "[#6:5]", S=f"{fg.X}")],
                    'pKa' : [None],
                    'Description' : 'Transmetalation', 
                },
                3 : {
                    'Templates' : [steps['Transmetalation']( "Ni", "Zn", f"[{fg.X}:4]", "[#6:5]")],
                    'pKa' : [None],
                    'Description' : '2nd Transmetalation', 
                },
                4 : {
                    'Templates' : [steps['Reductive elimination']['Normal']('Ni', '#6', '#6')],
                    'pKa' : [None],
                    'Description' : 'Reductive elimination to Ni(0)'

                },
                **core_cat_cycle([5, 6, 7], "Ni", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2;!$([Ni]-[Cl,Br,I,O])")
            }
        },

        # 'NiX2 and amidine' : {
        #     'Reagent': [f'[Ni](-[{fg.X}])(-[{fg.X}])',  fg.amidine_ligand , '[Zn;D0]',],
        #     'Exclude_reagent': [fg.amidine_ligand,],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Ni", fg.p_ligand , c = 0),
        #         2 : {
        #             'Templates' : [f'[Zn;D0:1].[Ni:2](-[{fg.X}:3])-[{fg.X}:4]>>[Zn:1](-[{fg.X}:3])-[{fg.X}:4].[Ni:2]'],
        #             'pKa' : [None],
        #             'Description' : 'Reduction of Ni(II)', 
        #         },
        #         3 : {
        #             'Templates' : [
        #                 '[#6:1]-[Br,I:2].[Zn;D0:3]>>[#6:1]-[Zn:3]-[Br,I:2]'
        #             ],
        #             'pKa' : [None],
        #             'Description' : 'In-situ organozinc preparation', 
        #         },
        #         **core_cat_cycle([4, 5, 6], "Ni", "Zn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2;!$([Ni]-[Cl,Br,I,O])")
        #     }
        # }
    },


    (
        'Bromo Stille reaction',
        'Chloro Stille reaction',
        'Iodo Stille reaction',
        'Triflyloxy Stille reaction',
        'Sulfoxy Stille reaction',
        # have to make sn partner more general
    ): 
    {
        **general_pd_coupling((lambda o, l: core_cat_cycle([o, o+1, o+2],"Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = l)), global_reagent=['[Sn]'], global_exclude=[fg.nhc, '[Fe;D0]']),
      
      

    },

    (
        'Bromo Heck reaction',
        'Chloro Heck reaction',
        'Iodo Heck reaction',
        'Triflyloxy Heck reaction',
        'Palladium arylhalide alkylation'
    ): {
        **general_pd_coupling(catalytic_cycles['Heck'], offset=0, global_reagent=[], global_exclude=[fg.nhc,'[Fe;D0]','[Zn,Cu]']),
    },

    (   
        'Goldberg reaction',
        'Bromo Goldberg reaction',
        'Chloro Goldberg reaction',
        'Iodo Goldberg reaction',
        'Fluoro Goldberg reaction',
    ): 
    {
        **general_cu_coupling(catalytic_cycles['Ullman C-N'], offset=0, partner="#7;D3", global_reagent=[f'[Cu;D1]-[{fg.X}]']),
        **general_cu_coupling(lambda o, r, l: {
            0: {
                'Templates' : [
                "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
                ],
                'pKa' : [None],
                'Description' : 'CuI2 Decomposition', 
            },
            **catalytic_cycles['Ullman C-N'](o, r, l)},
            offset = 1, partner = "#7;D3", modifier="CuI2", global_reagent = [f'[I]-[Cu;D2]-[I]']
        ),
        **general_cu_coupling(lambda o, r, l: {
            0: {
                'Templates' : [
                    '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                ],
                'pKa' : [None],
                'Description' : 'Tautomerization'
            },
            **catalytic_cycles['Ullman C-N'](o, r, l)},
            offset = 1, partner = "#7;D3", modifier="Tautomerization", global_reagent = [f'[Cu;D1]-[{fg.X}]']),
        **general_cu_coupling(lambda o, r, l: { 
            0 : {
                'Templates' : [
                "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
                ],
                'pKa' : [None],
                'Description' : 'CuI2 Decomposition', 
            },
            1: {
                'Templates' : [
                    '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                ],
                'pKa' : [None],
                'Description' : 'Tautomerization'
            },
            **catalytic_cycles['Ullman C-N'](o, r, l)},
            offset = 2, partner = "#7;D3", modifier="Tautomerization CuI2", global_reagent = [f'[I]-[Cu;D2]-[I]'])
    },

    (
        'Ullmann condensation'
    ):
    {
        **general_cu_coupling(catalytic_cycles['Ullman C-O'], offset=0, partner="#8", global_reagent=[f'[Cu;D1]-[{fg.X}]']),
        **general_cu_coupling(lambda o, r, l: {
            0: {
                'Templates' : [
                "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
                ],
                'pKa' : [None],
                'Description' : 'CuI2 Decomposition', 
            },
            **catalytic_cycles['Ullman C-O'](o, r, l)},
            offset = 1, partner = "#8", modifier="CuI2", global_reagent = [f'[I]-[Cu;D2]-[I]']
        ),
        **general_cu_coupling(lambda o, r, l: {
            0: {
                'Templates' : [
                    '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                ],
                'pKa' : [None],
                'Description' : 'Tautomerization'
            },
            **catalytic_cycles['Ullman C-O'](o, r, l)},
            offset = 1, partner = "#8", modifier="Tautomerization", global_reagent = [f'[Cu;D1]-[{fg.X}]']),
        **general_cu_coupling(lambda o, r, l: { 
            0 : {
                'Templates' : [
                "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
                ],
                'pKa' : [None],
                'Description' : 'CuI2 Decomposition', 
            },
            1: {
                'Templates' : [
                    '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                ],
                'pKa' : [None],
                'Description' : 'Tautomerization'
            },
            **catalytic_cycles['Ullman C-O'](o, r, l)},
            offset = 2, partner = "#8", modifier="Tautomerization CuI2", global_reagent = [f'[I]-[Cu;D2]-[I]']),
        # 'Cu(I)-X ligandless':{
        #     'Reagent': [f'[Cu]-[{fg.X}]'],
        #     'Exclude_reagent': [f'[{fg.n_n_ligand}]',f'[{fg.n_o_ligand}]','[O]=[C]-[C]-[C]=[O]','[Pd]'],
        #     'Stages': {
        #         0: {
        #         'Templates' : [
        #             "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
        #         ],
        #         'pKa' : [None],
        #         'Description' : 'CuI2 Decomposition', 
        #         },
        #         **catalytic_cycles['Ullman C-O'](1, "#8", 1)
        #     }
        # },
        # 'Cu(I)-X with ligand':{
        #     'Reagent': [f'[Cu]-[{fg.X}]',f'[{fg.n_n_ligand}]'],
        #     'Exclude_reagent': ['[Pd]'],
        #     'Stages': {
        #         0: {
        #         'Templates' : [
        #             "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
        #         ],
        #         'pKa' : [None],
        #         'Description' : 'CuI2 Decomposition', 
        #         },
        #         **catalyst_preparation['Attach ligands'](1, "Cu", f"#7;{fg.n_n_ligand}", c = 1),
        #         **catalytic_cycles['Ullman C-O'](3, "#8", 3)}
        # },
        # 'Cu(I)-X with n,o ligand':{
        #     'Reagent': [f'[Cu]-[{fg.X}]',f'[{fg.n_o_ligand}]'],
        #     'Exclude_reagent': ['[Pd]'],
        #     'Stages': {
        #         0: {
        #         'Templates' : [
        #             "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
        #         ],
        #         'pKa' : [None],
        #         'Description' : 'CuI2 Decomposition', 
        #         },
        #         1 : {
        #             'Templates': ["[O;H1:1]-[C:2](=[O:3])-[C:4]-[N:5]>>[O;-1;H0:1]-[C:2](=[O:3])-[C:4]-[N:5]",
        #                           "[O;H1:1]-[c:2]:[c:4]:[n:5]>>[O;H0;-1:1]-[c:2]:[c:4]:[n:5]",
        #                           "[O;H1:1]-[C:2](=[O:3])-[c:4]:[n:5]>>[O;H0;-1:1]-[C:2](=[O:3])-[c:4]:[n:5]"
        #                           ],
        #             'pKa': [{'B' : 3}, {'B' : 3}, {'B' : 3}],
        #             'Description': 'Deprotonate'
        #         },
        #         **catalyst_preparation['Attach ligands'](2, "Cu", f"{fg.o_n_ligand}", c = 1, lig_2=f"{fg.n_o_ligand}"),
        #         **catalytic_cycles['Ullman C-O'](4, "#8", 3)}
        # },
        # 'Cu(I)-X with diketone':{
        #     'Reagent': [f'[Cu]-[{fg.X}]','[O]=[C]-[C]-[C]=[O]'],
        #     'Exclude_reagent': ['[Pd]'],
        #     'Stages': {
        #         0: {
        #         'Templates' : [
        #             "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
        #         ],
        #         'pKa' : [None],
        #         'Description' : 'CuI2 Decomposition', 
        #         },
        #         1 : {
        #             'Templates': ["[O:1]=[C:2]-[C;H1:3]-[C:4]=[O:5]>>[O;-1:1]-[C:2]=[C;H0:3]-[C:4]=[O:5]",
        #                           "[O:1]=[C:2]-[C;H2:3]-[C:4]=[O:5]>>[O;-1:1]-[C:2]=[C;H1:3]-[C:4]=[O:5]"
        #                           ],
        #             'pKa': [{'B' : 3},{'B' : 3}],
        #             'Description': 'Deprotonate'
        #         },
        #         **catalyst_preparation['Attach ligands'](2, "Cu", f"O;-1;$([O]-[C]=[C]-[C]=[O])", c = 1, lig_2=f"O;$([O]=[C]-[C]=[C]-[O]-[Cu])"),
        #         **catalytic_cycles['Ullman C-O'](4, "#8", 3)}
        # },
    },

    ('Alkyne to alkane hydrogenation', 
    ) : 
    {
        'Pd and H2' : {
            'Reagent': ['[Pd;D0;+0]','[H]-[H]'],
            'Exclude_reagent': None,
            'Stages': catalytic_cycles['Alkyne hydrogenation'](0, "Pd")
            
        }, 
    },
    
    ('Alkene hydrogenation', 
    ) : 
    {
        'Pd and H2' : {
            'Reagent': ['[Pd;D0;+0]','[H]-[H]'],
            'Exclude_reagent': None,
            'Stages': catalytic_cycles['Alkene hydrogenation'](0, "Pd")
            
        }, 
    },

    ('Bromo Buchwald-Hartwig amination',
    'Chloro Buchwald-Hartwig amination',
    'Iodo Buchwald-Hartwig amination',
    'Triflyloxy Buchwald-Hartwig amination',) : 
    {
        **general_pd_coupling(catalytic_cycles['Buchwald-Hartwig amination'], offset=0, global_reagent=[], global_exclude=[fg.nhc, '[Fe;D0]'])
        # 'Pd(0) with no ligand' : {
        #   'Reagent': ['[Pd;D0;+0]'],
        #     'Exclude_reagent': ['[Fe;D0]', f'[{fg.p_ligand}]', fg.nhc],
        #     'Stages': {
        #         **catalytic_cycles['Buchwald-Hartwig amination'](0, n_l = 0, l = ";D0")},
        # },
        # 'Pd(0) with ligand' : {
        #   'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]'],
        #     'Exclude_reagent': ['[Fe;D0]',fg.nhc],
        #     'Stages': {
        #         0 : {
        #             'Templates': [steps["Ligand association"]["First"]("Pd", fg.p_ligand, 0, "!$([P]-[Pd])"),],
        #             'pKa': [None],
        #             'Description': 'Attach first ligand'
        #         },
        #         **catalytic_cycles['Buchwald-Hartwig amination'](1)},
        # },
        # 'Pd(0) with two, three, or four ligands' : {
        #     'Reagent': [f'[Pd;D2,D3,D4;$([Pd](-[*])-[*]),$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'],
        #     'Exclude_reagent': ['[Fe;D0]',fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Remove three ligands'](0, "Pd", fg.p_ligand),
        #         **catalytic_cycles['Buchwald-Hartwig amination'](3)},
        # },
        # 'Pd Buchwald precatalyst' : {
        #     'Reagent': [f'[{fg.anion}]-[Pd]-[#6]:[#6]-[#6]:[#6]-[#7;H2,H1]'],
        #     'Exclude_reagent': ['[Fe;D0]',fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Buchwald precatalyst'](0),
        #         **catalytic_cycles['Buchwald-Hartwig amination'](2)},
        # },
        # 'Pd(OAc)2 monophosphine' : {
        #     'Reagent': [fg.pd_oac_2, f'[{fg.p_ligand}]'],
        #     'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine,fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
        #         **catalyst_preparation['Pd(OAc)2 without water']['Monophosphine'](2, fg.p_ligand),
        #         7: {
        #             'Templates': [steps["Ligand dissociation"]["Third"]("Pd", fg.p_ligand),],
        #             'pKa': [None],
        #             'Description': 'Dissociate ligand'
        #         },
        #         **catalytic_cycles['Buchwald-Hartwig amination'](8)
        #     },
        # },
        # 'Pd(OAc)2 bisphosphine' : {
        #     'Reagent': [fg.pd_oac_2,fg.bisphosphine],
        #     'Exclude_reagent': ['[Fe;D0]',fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
        #         **catalyst_preparation['Pd(OAc)2 without water']['Bisphosphine'](2, fg.p_ligand),
        #         **catalytic_cycles['Buchwald-Hartwig amination'](6)
        #     },
        # }
    },

    ('Bromo Miyaura boration',
    'Chloro Miyaura boration',
    'Iodo Miyaura boration',
    'Triflyloxy Miyaura boration',) : {
        # 'Pd(0) with no ligand' : {
        #   'Reagent': ['[Pd;D0;+0]'],
        #     'Exclude_reagent': ['[Fe;D0]', f'[{fg.p_ligand}]'],
        #     'Stages': {
        #         **catalytic_cycles['Miyaura borylation'](0, n_l = 0, l = ";D0")},
        # },
        # 'Pd(0) with ligand' : {
        #   'Reagent': ['[Pd;D0;+0]',f'[{fg.p_ligand}]'],
        #     'Exclude_reagent': ['[Fe;D0]'],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 0),
        #         **catalytic_cycles['Miyaura borylation'](2)},
        # },
        **general_pd_coupling(catalytic_cycles['Miyaura borylation'], offset=0, global_reagent=[], global_exclude=[fg.nhc, '[Fe;D0]'])
    },

    # ('Simmons-Smith-Furukawa reaction') : {
    #     'Furukawa' : {
    #         'Reagent': ['[C;H3]-[C;H2]-[Zn;D2;+0]-[C;H2]-[C;H3]',f'[C;H2](-[I,Cl])-[I]'],
    #         'Exclude_reagent' : ['[C;$([C](-[F])(-[F])-[F])]-[C](=[O])-[O;H1]'],
    #         'Stages': {
    #             0 : {
    #                 'Templates': ['[Zn;D2;+0:1]-[#6:2].[C;H2:3](-[I,Cl:5])-[I:4]>>[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:5].[#6:2]-[I:4]'
    #                 ],
    #             'pKa': [None],
    #             'Description': 'Zinc-halogen exchange',
    #             },
    #             1 : {
    #                 'Templates': ['[#6:2]-[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:4].[#6:5]=[#6:6]>>[I,Cl:4]-[Zn:1]-[#6:2].[#6:5]1-[#6:6]-[C;H2:3]-1',
    #                             '[#6:2]-[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:4].[#6:5]:[#6:6]>>[I,Cl:4]-[Zn:1]-[#6:2].[#6:5]1-[#6:6]-[C;H2:3]-1'
    #                 ],
    #             'pKa': [None],
    #             'Description': 'Cyclopropanation',
    #             },
    #         }

    #     },
    #     'Shi' : {
    #         'Reagent': ['[C;H3]-[C;H2]-[Zn;D2;+0]-[C;H2]-[C;H3]',f'[C;H2](-[I,Cl])-[I]','[C;$([C](-[F])(-[F])-[F])]-[C](=[O])-[O;H1]'],
    #         'Exclude_reagent' : None,
    #         'Stages': {
    #             0 : {
    #                 'Templates': ['[Zn;D2;+0:1]-[#6;H2:2].[C;$([C](-[F])(-[F])-[F]):6]-[C:5](=[O:4])-[O;H1:3]>>[#6;H3:2].[C;$([C](-[F])(-[F])-[F]):6]-[C:5](=[O:4])-[O;H0:3]-[Zn;D2;+0:1]'
    #                 ],
    #             'pKa': [None],
    #             'Description': 'TFA for ethyl',
    #             },
    #             1 : {
    #                 'Templates': ['[O:6]-[Zn;D2;+0:1]-[#6:2].[C;H2:3](-[I,Cl:5])-[I:4]>>[O:6]-[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:5].[#6:2]-[I:4]'
    #                 ],
    #             'pKa': [None],
    #             'Description': 'Zinc-halogen exchange',
    #             },
    #             2 : {
    #                 'Templates': ['[O:2]-[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:4].[#6:5]=[#6:6]>>[I,Cl:4]-[Zn:1]-[O:2].[#6:5]1-[#6:6]-[C;H2:3]-1',
    #                               '[O:2]-[Zn;D2;+0:1]-[C;H2:3]-[I,Cl:4].[#6:5]:[#6:6]>>[I,Cl:4]-[Zn:1]-[O:2].[#6:5]1-[#6:6]-[C;H2:3]-1'
    #                 ],
    #             'pKa': [None],
    #             'Description': 'Cyclopropanation',
    #             },
    #         },
            
    #     }
    # },

    ('Nitro to amino') : {
        'Pd and H2' : {
            'Reagent': ['[Pd;D0;+0]','[H]-[H]'],
            'Exclude_reagent': None,
            'Stages': catalytic_cycles["Nitro to amine"](0,"Pd")
        }, 
    },

    ('Azido to amino') : 
    {
        'Pd and H2' : {
            'Reagent': ['[Pd;D0;+0]','[H]-[H]'],
            'Exclude_reagent': None,
            'Stages': catalytic_cycles["Azide to amine"](0,"Pd")
            
        }, 
    },

    ('N-Bn deprotection','O-Bn deprotection','CO2H-Bn deprotection', 'N-Cbz deprotection') : 
    
    {
        'Pd and H2' : {
            'Reagent' : ['[Pd;D0;+0]','[H]-[H]'],
            'Exclude_reagent' : ['[N]-[C](=[O])-[O]-[C]-c1ccccc1'],
            'Stages': catalytic_cycles['Debenzylation'](0)
        },
        'Pd and H2 and Cbz' : {
            'Reagent' : ['[Pd;D0;+0]','[H]-[H]','[N]-[C](=[O])-[O]-[C]-c1ccccc1'],
            'Exclude_reagent' : None,
            'Stages': {
                **catalytic_cycles['Debenzylation'](0),
                3 : {
                    'Templates': ['[N:1]-[C:2](=[O:3])-[O:4]-[H:5]>>[N:1]-[H:5].[O:3]=[C:2]=[O:4]'
                    ],
                    'pKa': [None],
                    'Description': 'Decarboxylation',
                }
            }
        }
    },

    ('Palladium-catalyzed direct arylation') :
    {
        **general_pd_coupling(catalytic_cycles['CMD C-H arylation'], offset=0, global_reagent=[], global_exclude=[fg.nhc,'[Fe;D0]','[Zn,Cu]']),
        'Pd(OAc)2 monophosphine ' : {
            'Reagent': [fg.pd_oac_2, f'[{fg.p_ligand}]'],
            'Exclude_reagent': ['[Fe;D0]', fg.bisphosphine, fg.nhc, "[Zn,Cu]"],
            'Stages': {
                **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
                **catalyst_preparation['Pd(OAc)2 without water']['Monophosphine'](2, fg.p_ligand),
                **catalytic_cycles['CMD C-H arylation'](7, ";D2;!$([Pd]-[Cl,Br,O])")
            },
        },
        **general_pd_coupling(catalytic_cycles['Heck'], offset=0, global_reagent=[], global_exclude=[fg.nhc,'[Fe;D0]','[Zn,Cu]'], modifier="Heck"),
       
    },

    ('Migita thioether synthesis') : {
        **general_pd_coupling(catalytic_cycles['Migita thioether'], offset = 0, global_reagent=['[S]'], global_exclude=[fg.nhc,'[Fe;D0]'])
    },

    ('Bromo alkoxycarbonylation', "Chloro alkoxycarbonylation", "Iodo alkoxycarbonylation") : {
        **general_pd_coupling(catalytic_cycles['Alkoxycarbonylation'], offset=0, global_reagent=['[C]#[O]'], global_exclude=[fg.nhc,'[Fe;D0]'])
    },
    ('Bromo aminocarbonylation', "Chloro aminocarbonylation", "Iodo aminocarbonylation") : {
        **general_pd_coupling(catalytic_cycles['Aminocarbonylation'], offset=0, global_reagent=['[C]#[O]'], global_exclude=[fg.nhc,'[Fe;D0]'])
    },

    ('Chan-Lam alkylamine coupling', 'Chan-Lam arylamine coupling') : {
        'Without oxygen' : {
            'Reagent': ['[Cu]'],
            'Exclude_reagent' : ['[O]=[O]'],
            'Stages': {
                **catalytic_cycles['Chan-Lam C-N'](0),
            }
        },
        'Without oxygen - taut' : {
            'Reagent': ['[Cu]'],
            'Exclude_reagent' : ['[O]=[O]'],
            'Stages': {
                0 : {
                    'Templates' : [
                        '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                    ],
                    'pKa' : [None],
                    'Description' : 'Tautomerization'
                },
                **catalytic_cycles['Chan-Lam C-N'](1),
            }
        },
        # 'With oxygen' : {
        #     'Reagent': ['[Cu]','[O]=[O]'],
        #     'Exclude_reagent' : None,
        #     'Stages': {
        #         **catalytic_cycles['Chan-Lam C-N'](0),
        #         3 : {

        #         }
        #     }
        # }
    },

    ('Chan-Lam ether coupling') : {
        'Without oxygen' : {
            'Reagent': ['[Cu]'],
            'Exclude_reagent' : ['[O]=[O]'],
            'Stages': {
                **catalytic_cycles['Chan-Lam C-O'](0),
            }
        },
        'Without oxygen - taut' : {
            'Reagent': ['[Cu]'],
            'Exclude_reagent' : ['[O]=[O]'],
            'Stages': {
                0 : {
                    'Templates' : [
                        '[n;H1:1]:[n;H0:2]>>[n;H0:1]:[n;H1:2]',
                    ],
                    'pKa' : [None],
                    'Description' : 'Tautomerization'
                },
                **catalytic_cycles['Chan-Lam C-O'](1),
            }
        },
    },

    ('Rosenmund von Braun cyanation') : {
        'CuCN' : {
            'Reagent' : ['[Cu]-[C]#[N]'],
            'Exclude_reagent': None,
            'Stages' : {
                **catalytic_cycles['Rosenmund'](0)
            }
        }
    },

    ('Dechlorination', 'Debromination', 'Deiodination') : {
        **general_pd_coupling(catalytic_cycles['Pd alcohol dehalogenation'], offset=0, global_reagent=['[O;H1;!$([O]-[C]=[O])]-[#6]'], global_exclude=['[C;H1](=[O])-[O]']),
        **general_pd_coupling(catalytic_cycles["Pd dehalogenation with hydrogen"], offset=0, global_reagent=['[H]-[H]'], global_exclude=['[Si;H1,H2,H3]','[C;H1](=[O])-[O]'], modifier=" hydrogen"),
        **general_pd_coupling(catalytic_cycles["Pd dehalogenation with hydrogen"], offset=0, global_reagent=['[Si;H1,H2,H3]'], global_exclude=['[C;H1](=[O])-[O]'], modifier=" silane"),
        **general_pd_coupling(catalytic_cycles["Pd formate hydrogenation"], offset=0, global_reagent=['[C;H1](=[O])-[O]'], global_exclude=['[Si;H1,H2,H3]'], modifier=" formate")
    },

    ('Hiyama coupling') : {
        **general_pd_coupling(catalytic_cycles['Hiyama'], offset=0, global_reagent=['[F]'], global_exclude=[fg.nhc, '[Fe;D0]', '[Cu]']),
        **general_pd_coupling(catalytic_cycles['Hiyama with copper'], offset=0, global_reagent=['[Cu]', '[F]'], global_exclude=[fg.nhc, '[Fe;D0]'], modifier="copper"),
        **general_pd_coupling(catalytic_cycles['Hiyama base activation with Cu'], offset=0, global_reagent=['[Cu]'], global_exclude=[fg.nhc, '[Fe;D0]', '[F]'], modifier="copper and base"),
        **general_pd_coupling(catalytic_cycles['Hiyama base activation'], offset=0, global_reagent=[], global_exclude=[fg.nhc, '[Fe;D0]', '[F]', '[Cu]'], modifier="base"),
         **general_pd_coupling(catalytic_cycles['Hiyama with copper, electron rich Si'], offset=0, global_reagent=['[Cu]'], global_exclude=[fg.nhc, '[Fe;D0]', '[F]'], modifier="electron rich silane"),
    },

    ('Keto alpha-arylation') : {
        **general_pd_coupling(catalytic_cycles['Alpha arylation'], offset=0, global_reagent=[], global_exclude=[fg.nhc, '[Fe;D0]'])
    },

    ('O-allyl deprotection') : {
        "With amine" : {
            'Reagent' : ['[N;H1,H2,H3]', f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'], 
            'Exclude_reagent' : None, 
            'Stages' : {
                **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
                2 : {
                    'Templates' : ['[Pd;D2;!$([Pd]-[Cl,Br,I,O]):1].[O:2]-[C:3]-[C:4]=[C:5]>>[O;-1:2].[Pd;+1:1]-[C:3]-[C:4]=[C:5]'],
                    'pKa' : [None],
                    'Description' : 'Oxidative addition'
                },
                3 : {
                    'Templates' : ['[N:1]-[H:2].[O;-1:3].[Pd;+1:4]-[C:5]-[C:6]=[C:7]>>[O;+0:3]-[H:2].[N:1]-[C:7]-[C:6]=[C:5].[Pd;+0:4]'],
                    'pKa' : [None],
                    'Description' : 'Amine attack'
                }
            }
        },
        "With alcohol" : {
            'Reagent' : ['[C;H3]-[O;H1]', f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'], 
            'Exclude_reagent' : None, 
            'Stages' : {
                **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
                2 : {
                    'Templates' : ['[Pd;D2;!$([Pd]-[Cl,Br,I,O]):1].[O:2]-[C:3]-[C:4]=[C:5]>>[O;-1:2].[Pd;+1:1]-[C:3]-[C:4]=[C:5]'],
                    'pKa' : [None],
                    'Description' : 'Oxidative addition'
                },
                3 : {
                    'Templates' : ['[C;H3:1]-[O;H1:2]>>[O;-1:2]-[C;H3:1]'],
                    'pKa' : [{'B' : 3}],
                    'Description' : 'Deprotonation'
                },
                4 : {
                    'Templates' : ['[O;-1:3].[Pd;+1:4]-[C:5]-[C:6]=[C:7]>>[O;+0:3]-[C:7]-[C:6]=[C:5].[Pd;+0:4]'],
                    'pKa' : [None],
                    'Description' : 'Alkoxide attack'
                },
                5 : {
                    'Templates' : ['[O;-1:3]>>[O;+0;H1:3]'],
                    'pKa' : [{'A' : 15}],
                    'Description' : 'Protonation of product'
                }
            }
        }
    }


    # # for testing:
    # ('Test') : {
    #     'Water' : {
    #         'Reagent' : ['[O;H2]'],
    #         'Exclude_reagent': None,
    #         'Stages': {
    #             0 : {
    #                 'Templates': ['[O;H2;+0:1].[O;H2;+0:2]>>[O;H1;-1:1].[O;H3;+1:2]'
    #                 ],
    #                 'pKa': [None],
    #                 'Description': 'Reductive elimination',
    #             },
    #             # 1 : {
    #             #     'Templates': ['[O;H2;+0:1].[O;H2;+0:2]>>[O;H1;-1:1].[O;H3;+1:2]'
    #             #     ],
    #             #     'pKa': [None],
    #             #     'Description': 'Reductive elimination',
    #             # }
    #         }
    #     },
    #     'Copper' : {
    #         'Reagent': [f'[I]-[Cu;D2]-[I]'],
    #             'Exclude_reagent': [fg.n_n_ligand],
    #             'Stages': {
    #                 0: {
    #                     'Templates' : [
    #                         "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D2:5]-[I:6]>>[Cu;D1:2]-[I:3].[Cu;D1:5]-[I:4].[I:1]-[I:6]"
    #                     ],
    #                     'pKa' : [None],
    #                     'Description' : 'CuI2 Decomposition', 
    #                 },
    #                 1: {
    #                     'Templates' : [
    #                         "[I:1]-[Cu;D2:2]-[I:3].[I:4]-[Cu;D1:5]>>[Cu;+2:2].[Cu;+1:5].[I;-1:1].[I;-1:3].[I;-1:4]"
    #                     ],
    #                     'pKa' : [None],
    #                     'Description' : 'Cu2I3 Formation', 
    #                 },
    #         }
    #     }
    # }


}


##### graveyard

  # 'Pd(0) ligandless': {
        #     'Reagent': ['[Pd;D0;+0]'],
        #     'Exclude_reagent': ['[P]'],
        #     'Stages': {
        #         **core_cat_cycle([0, 1, 2], "Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D0")
        #     },
        # },
        # 'Pd(0) with ligand': {
        #     'Reagent': ['[Pd;D0;+0]', f'[{fg.p_ligand}]'],
        #     'Exclude_reagent': ['[Fe;D0]', fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])"),
        #         **core_cat_cycle([2, 3, 4], "Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2")
        #     }
        # },
        # 'Pd with three or four ligands': {
        #     'Reagent': [f'[Pd;D3,D4;$([Pd](-[*])(-[*])-[*]),$([Pd](-[*])(-[*])(-[*])-[*]);!$([Pd][{fg.anion}]);+0]'],
        #     'Exclude_reagent': ['[Fe;D0]', fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Remove ligands'](0, "Pd", fg.p_ligand),
        #         **core_cat_cycle([2, 3, 4], "Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D2")
        #     },
        # },
        # 'Pd(OAc)2 monophosphine' : {
        #     'Reagent': ['[C;H3]-[C](=[O])-[O]-[Pd;+0;D2]-[O]-[C](=[O])-[C;H3]', f'[{fg.p_ligand}]'],
        #     'Exclude_reagent': [fg.nhc, '[Fe;D0]', fg.bisphosphine],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
        #         **catalyst_preparation['Pd(OAc)2 without water']['Monophosphine'](2, fg.p_ligand),
        #         7: {
        #             'Templates': [steps["Ligand dissociation"]["Third"]("Pd", fg.p_ligand),],
        #             'pKa': [None],
        #             'Description': 'Dissociate ligand'
        #         },
        #         **core_cat_cycle([8, 9, 10], "Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D1")
        #     }
        # },
        # 'Pd(OAc)2 bisphosphine' : {
        #     'Reagent': ['[C;H3]-[C](=[O])-[O]-[Pd;+0;D2]-[O]-[C](=[O])-[C;H3]', fg.bisphosphine],
        #     'Exclude_reagent': ['[Fe;D0]',fg.nhc],
        #     'Stages': {
        #         **catalyst_preparation['Attach ligands'](0, "Pd", fg.p_ligand, "!$([P]-[Pd])", c = 2),
        #         **catalyst_preparation['Pd(OAc)2 without water']['Bisphosphine'](2, fg.p_ligand),
        #         **core_cat_cycle([6, 7, 8], "Pd", "Sn", f"[{fg.X}:4]", "[#6:5]", "#6", "#6", l1 = ";D1")
        #     },
        # }