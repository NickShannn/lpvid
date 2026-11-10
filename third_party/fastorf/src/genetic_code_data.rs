// 由 tools/gen_gc_tables.py 从 NCBI gc.prt (v4.6) 生成, 勿手改。
// 密码子索引 = b1*16 + b2*4 + b3, 碱基编码 A=0,C=1,G=2,T=3。
pub struct CodeTable {
    pub id: u32,
    pub name: &'static str,
    pub aa: &'static str,     // 64 个 IUPAC 氨基酸, * = stop
    pub starts: &'static str, // 64 位标志, 'M' = initiator
}

pub static TABLES: &[CodeTable] = &[
    CodeTable { id: 1, name: "Standard", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", starts: "--------------M---------------M-----------------*-*-----*-----M-" },
    CodeTable { id: 2, name: "Vertebrate Mitochondrial", aa: "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "--------*-*-MMMM------------------------------M-*-*-------------" },
    CodeTable { id: 3, name: "Yeast Mitochondrial", aa: "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "------------M-M-------------------------------M-*-*-------------" },
    CodeTable { id: 4, name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate
 Mitochondrial; Mycoplasma; Spiroplasma", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "------------MMMM--------------M---------------M-*-*---------M-M-" },
    CodeTable { id: 5, name: "Invertebrate Mitochondrial", aa: "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "------------MMMM------------------------------M-*-*-----------M-" },
    CodeTable { id: 6, name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF", starts: "--------------M-----------------------------------------*-------" },
    CodeTable { id: 9, name: "Echinoderm Mitochondrial; Flatworm Mitochondrial", aa: "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "--------------M-------------------------------M-*-*-------------" },
    CodeTable { id: 10, name: "Euplotid Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF", starts: "--------------M---------------------------------*-*-------------" },
    CodeTable { id: 11, name: "Bacterial, Archaeal and Plant Plastid", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", starts: "------------MMMM--------------M---------------M-*-*-----*-----M-" },
    CodeTable { id: 12, name: "Alternative Yeast Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", starts: "--------------M---------------M-----------------*-*-----*-------" },
    CodeTable { id: 13, name: "Ascidian Mitochondrial", aa: "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "------------M-M-------------------------------M-*-*-----------M-" },
    CodeTable { id: 14, name: "Alternative Flatworm Mitochondrial", aa: "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF", starts: "--------------M-----------------------------------*-------------" },
    CodeTable { id: 15, name: "Blepharisma Macronuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF", starts: "--------------M---------------------------------*-------*-------" },
    CodeTable { id: 16, name: "Chlorophycean Mitochondrial", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF", starts: "--------------M---------------------------------*-------*-------" },
    CodeTable { id: 21, name: "Trematode Mitochondrial", aa: "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "--------------M-------------------------------M-*-*-------------" },
    CodeTable { id: 22, name: "Scenedesmus obliquus Mitochondrial", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF", starts: "--------------M---------------------------------*---*---*-------" },
    CodeTable { id: 23, name: "Thraustochytrium Mitochondrial", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF", starts: "--------------MM------------------------------M-*-*-----*---*---" },
    CodeTable { id: 24, name: "Rhabdopleuridae Mitochondrial", aa: "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF", starts: "--------------M---------------M---------------M-*-*-----------M-" },
    CodeTable { id: 25, name: "Candidate Division SR1 and Gracilibacteria", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSGCWCLFLF", starts: "--------------M-------------------------------M-*-*-----------M-" },
    CodeTable { id: 26, name: "Pachysolen tannophilus Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLALEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF", starts: "--------------M---------------M-----------------*-*-----*-------" },
    CodeTable { id: 27, name: "Karyorelict Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF", starts: "--------------M-----------------------------------------*-------" },
    CodeTable { id: 28, name: "Condylostoma Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSSWCWCLFLF", starts: "--------------M---------------------------------*-*-----*-------" },
    CodeTable { id: 29, name: "Mesodinium Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYYYYSSSS*CWCLFLF", starts: "--------------M-----------------------------------------*-------" },
    CodeTable { id: 30, name: "Peritrich Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSS*CWCLFLF", starts: "--------------M-----------------------------------------*-------" },
    CodeTable { id: 31, name: "Blastocrithidia Nuclear", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVEYEYSSSSWCWCLFLF", starts: "--------------M---------------------------------*-*-------------" },
    CodeTable { id: 32, name: "Balanophoraceae Plastid", aa: "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YWYSSSS*CWCLFLF", starts: "------------MMMM--------------M---------------M-*-------*-----M-" },
    CodeTable { id: 33, name: "Cephalodiscidae Mitochondrial", aa: "KNKNTTTTSSKSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF", starts: "--------------M---------------M---------------M---*-----------M-" },
];
