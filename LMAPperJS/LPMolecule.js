	//From Lmapper, https://sourceforge.net/p/spm-and-mol-viewer/code/HEAD/tree/LP_SPM_MolView/LPAtom.vb
const PeriodicTableNames = ["H", "He",
	 "Li", "Be", "B", "C", "N", "O", "F", "Ne",
	 "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
	 "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
	 "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
	 "Cs", "Ba",
		"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
				"Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
	 "Fr", "Ra",
		"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
				"Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
	];

const PeriodicTableCovalentRadius_Angst=
	[0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03, 1.76, 1.7, 1.6, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.2, 1.19, 1.2, 1.2, 1.16, 2.2, 1.95, 1.9, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44, 1.42, 1.39, 1.39, 1.38, 1.39, 1.4, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.9, 1.87, 1.87, 1.75, 1.7, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32, 1.45, 1.46, 1.48, 1.4, 1.5, 1.5, 2.6, 2.21, 2.15, 2.06, 2, 1.96, 1.9, 1.87, 1.8, 1.69, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6];

const PeriodicTableVDWRadius_Angst= [1.1, 1.4, 1.81, 1.53, 1.92, 1.7, 1.55, 1.52, 1.47, 1.54, 2.27, 1.73, 1.84, 2.1, 1.8, 1.8, 1.75, 1.88, 2.75, 2.31, 2.3, 2.15, 2.05, 2.05, 2.05, 2.05, 2, 2, 2, 2.1, 1.87, 2.11, 1.85, 1.9, 1.83, 2.02, 3.03, 2.49, 2.4, 2.3, 2.15, 2.1, 2.05, 2.05, 2, 2.05, 2.1, 2.2, 2.2, 1.93, 2.17, 2.06, 1.98, 2.16, 3.43, 2.68, 2.5, 2.48, 2.47, 2.45, 2.43, 2.42, 2.4, 2.38, 2.37, 2.35, 2.33, 2.32, 2.3, 2.28, 2.27, 2.25, 2.2, 2.1, 2.05, 2, 2, 2.05, 2.1, 2.05, 1.96, 2.02, 2.07, 1.97, 2.02, 2.2, 3.48, 2.83, 2, 2.4, 2, 2.3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];

//Colors in hex string format
const PeriodicTableColor=
	["FFFFFF",  //H
	 "D9FFFF",  //He
	 "CC80FF",  //Li
	 "C2FF00",  //Be
	 "FFB5B5",  //B
	 "909090",  //C
	 "3050F8",  //N
	 "FF0D0D",  //O
	 "90E050",  //F
	 "B3E3F5",  //Ne
	 "AB5CF2",  //Na
	 "8AFF00",  //Mg
	 "BFA6A6",  //Al
	 "F0C8A0",  //Si
	 "FF8000",  //P
	 "FFFF30",  //S
	 "1FF01F",  //Cl
	 "80D1E3",  //Ar
	 "8F40D4",  //K
	 "3DFF00",  //Ca
	 "E6E6E6",  //Sc
	 "BFC2C7",  //Ti
	 "A6A6AB",  //V
	 "8A99C7",  //Cr
	 "9C7AC7",  //Mn
	 "E06633",  //Fe
	 "F090A0",  //Co
	 "50D050",  //Ni
	 "C88033",  //Cu
	 "7D80B0",  //Zn
	 "C28F8F",  //Ga
	 "668F8F",  //Ge
	 "BD80E3",  //As
	 "FFA100",  //Se
	 "A62929",  //Br
	 "5CB8D1",  //Kr
	 "702EB0",  //Rb
	 "00FF00",  //Sr
	 "94FFFF",  //Y
	 "94E0E0",  //Zr
	 "73C2C9",  //Nb
	 "54B5B5",  //Mo
	 "3B9E9E",  //Tc
	 "248F8F",  //Ru
	 "0A7D8C",  //Rh
	 "006985",  //Pd
	 "C0C0C0",  //Ag
	 "FFD98F",  //Cd
	 "A67573",  //In
	 "668080",  //Sn
	 "9E63B5",  //Sb
	 "D47A00",  //Te
	 "940094",  //I
	 "429EB0",  //Xe
	 "57178F",
	 "00C900",
	 "70D4FF",
	 "FFFFC7",
	 "D9FFC7",
	 "C7FFC7",
	 "A3FFC7",
	 "8FFFC7",
	 "61FFC7",
	 "45FFC7",
	 "30FFC7",
	 "1FFFC7",
	 "00FF9C",
	 "00E675",
	 "00D452",
	 "00BF38",
	 "00AB24",
	 "4DC2FF",
	 "4DA6FF",
	 "2194D6",
	 "267DAB",
	 "266696",
	 "175487",
	 "D0D0E0",
	 "FFD123",
	 "B8B8D0",
	 "A6544D",
	 "575961",
	 "9E4FB5",
	 "AB5C00",
	 "754F45",
	 "428296",
	 "420066",
	 "007D00",
	 "70ABFA",
	 "00BAFF",
	 "00A1FF",
	 "008FFF",
	 "0080FF",
	 "006BFF",
	 "545CF2",
	 "785CE3",
	 "8A4FE3",
	 "A136D4",
	 "B31FD4",
	 "B31FBA",
	 "B30DA6",
	 "BD0D87",
	 "C70066",
	 "CC0059",
	 "D1004F",
	 "D90045",
	 "E00038",
	 "E6002E",
	 "EB0026"
	];


class LPAtom{
	constructor(elName, x_Angst , y_Angst , z_Angst ,  ser_number_string){
		//this.elementName= elName;

		this.elNumber=-1; //undefined

		this.x_Angst = x_Angst;
		this.y_Angst = y_Angst;
		this.z_Angst = z_Angst;
		this.ser_number_string=ser_number_string; //serial number or id
		
		let element_found_index = PeriodicTableNames.findIndex(
			function(pt_elname) { return elName.trim().toLowerCase()== pt_elname.toLowerCase() } 
		)
		if (element_found_index>=0){
			this.elNumber = element_found_index+1;
		}

	}

	getElementName(){
		if (this.elNumber < 1){
			return "";
		}else{
			return PeriodicTableNames[this.elNumber-1];
		}
	}

	getElementColorRGBHexString(){
		return PeriodicTableColor[this.elNumber-1];
	}

	getCovalentRadius_Angst(){
		return PeriodicTableCovalentRadius_Angst[this.elNumber-1];
	}
	getVDWRadius_Angst(){
		return PeriodicTableVDWRadius_Angst[this.elNumber-1];
	}
	/*
	getElementNumber(){
		//Make sure name is valid
		let b_findElNumber=false;

		if (this.elNumber==-1){
			b_findElNumber=true;
		}else{
			if (this.PeriodicTableNames[elNumber-1].toLowerCase != this.elementName ) b_findElNumber=true;
		}

		if (b_findElNumber){
			//Use the function Array.prototype.findIndex() , that uses a callback function
			let element_found_index = PeriodicTableNames.findIndex(
				function(pt_elname) { return this.elementName.trim().toLowerCase()== pt_elname.toLowerCase() } 
			)
			if (element_found_index>=0){
				this.elNumber = element_found_index+1;
			}
		}

		return this.elNumber;
	}
	*/
	clone(){
		let cl_atom=new LPAtom(
			this.getElementName(),
			this.x_Angst ,
			this.y_Angst,
			this.z_Angst,
			this.ser_number_string
		);
		return cl_atom;
	}
}

class LPBond{
	constructor( at1, at2, bondOrder=1){
		//Should check 
		this.atom1=at1;
		this.atom2=at2;
		this.bondOrder=bondOrder;
	}
	//Cannot have a clone method becasue it strongly depends on the atom objects
}

class LPMolecule{

	constructor(){
		this.listAtoms=[];
		this.listBonds=[];
		this.moleculeName="";
		this.filename="";
	}

	addAtom( lPAtomObject ){
		this.listAtoms.push(lPAtomObject);
	}

	addBond(lpBondObj){
		this.listBonds.push(lpBondObj);
	}

	static getTextFromArrayBuffer(arrbuff){
		let u8=new Uint8Array(arrbuff);
		let s= String.fromCharCode.apply( null, u8 );

		return s;
	}

	static getLPMoleculeFromArrBuffCMLFile(arrbuff0){
		//Creates a LPMolecule object by parsing a blob as
		//Chemical Markup Language (CML),
		//a type of XML

		
		var m_molecule=null;
		
		/*
		var enc = new TextDecoder("utf-8");
		var dec_s_ab = enc.decode(arrbuff0);
		*/

		var dec_s_ab= LPMolecule.getTextFromArrayBuffer(arrbuff0);
	
		let parser = new DOMParser();
		let doc_DOM = parser.parseFromString(dec_s_ab, "text/xml");

		//Gets the first molecule only
		var doc_DOM_molecules = doc_DOM.getElementsByTagName("molecule");

		if ( doc_DOM_molecules.length >0 ){
			let doc_molecule = doc_DOM_molecules[0];
			
			m_molecule=new LPMolecule();
			m_molecule.moleculeName= doc_molecule.getAttribute("id");

			let bHasCrystalParams=false;
			let crystal_a_angst=0;
			let crystal_b_angst=0;
			let crystal_c_angst=0;
			let crystal_alpha_rad=0;
			let crystal_beta_rad=0;
			let crystal_gamma_rad=0;

			let doc_mol_crystal = doc_molecule.getElementsByTagName("crystal");
			if (doc_mol_crystal.length > 0){
				bHasCrystalParams = true;

				//Parse crystal parameters
				
				let crystal_scalars = doc_mol_crystal[0].getElementsByTagName("scalar");

				for (let i=0 ; i< crystal_scalars.length ; i++){
					let c_scalar = crystal_scalars[i];

					//Get the value
					let v = parseFloat(c_scalar.innerHTML);

					//Correction to the value
					let corr=1.0;

					switch ( c_scalar.getAttribute( "units") ){
						case "units:angstrom":
							corr = 1.0;
							break;
						case "units:nm":
							corr = 10.0;
							break;
						case "units:degree":
							corr = Math.PI / 180;
							break;
						case ("units:rad" || "units:radians"):
							corr = 1.0;
							break;
					}

					//Apply correction npw
					v = v*corr;

					switch ( c_scalar.getAttribute("title") ){
						case "a":
							crystal_a_angst=v;
							break;
						case "b":
							crystal_b_angst=v;
							break;
						case "c":
							crystal_c_angst=v;
							break;
						case "alpha":
							crystal_alpha_rad=v;
							break;
						case "beta":
							crystal_beta_rad=v;
							break;
						case "gamma":
							crystal_gamma_rad=v;
							break;
					}

				}
			}
			let atoms=doc_molecule.getElementsByTagName("atomArray")[0].getElementsByTagName("atom"); //Gets list of atoms
			
			for (let i=0; i<atoms.length ; i++){
				let atom= atoms[i];
				let atom_id= atom.getAttribute("id"); //or sernumber
			
				let atom_type= atom.getAttribute("elementType");
				
				let atom_x3=0;
				let atom_y3=0;
				let atom_z3=0;

				if ( bHasCrystalParams ) {
					//if it is crystal these fractional values will be present
					let atom_xFract =  parseFloat(atom.getAttribute("xFract"));
					let atom_yFract =  parseFloat(atom.getAttribute("yFract"));
					let atom_zFract =  parseFloat(atom.getAttribute("zFract"));

					//Convert crystal fractional coordinates to Cartesian coordinates
					// https://en.wikipedia.org/wiki/Fractional_coordinates#Conversion_to_Cartesian_coordinates

					atom_x3 = crystal_a_angst * atom_xFract +
                        crystal_b_angst * Math.cos(crystal_gamma_rad) * atom_yFract +
						 crystal_c_angst * Math.cos( crystal_beta_rad ) * atom_zFract ;
						 
					atom_y3 = crystal_b_angst * Math.sin(crystal_gamma_rad) * atom_yFract +
						crystal_c_angst * ((Math.cos(crystal_alpha_rad) - Math.cos(crystal_beta_rad) * Math.cos(crystal_gamma_rad)) / Math.sin(crystal_gamma_rad)) * atom_zFract ;
					
                    let vol0 = crystal_a_angst * crystal_b_angst * crystal_c_angst * Math.sqrt(1 - Math.cos(crystal_alpha_rad) ** 2 -
                        Math.cos(crystal_beta_rad) ** 2 - Math.cos(crystal_gamma_rad) ** 2 +
                        2 * Math.cos(crystal_alpha_rad) * Math.cos(crystal_beta_rad) * Math.cos(crystal_gamma_rad)) ;
					
					atom_z3 = vol0 / crystal_a_angst / crystal_b_angst / Math.sin( crystal_gamma_rad ) * atom_zFract ;

				}else{
					atom_x3= parseFloat(atom.getAttribute("x3"));
					atom_y3= parseFloat(atom.getAttribute("y3"));
					atom_z3= parseFloat(atom.getAttribute("z3"));
				}

				let at0 = new LPAtom(atom_type , atom_x3, atom_y3, atom_z3, atom_id);
				//Add atom to molecule
				m_molecule.addAtom(at0);

			}

			let bonds=doc_molecule.getElementsByTagName("bondArray")[0].getElementsByTagName("bond"); //Gets list of bonds

			for (let i=0; i<bonds.length ; i++){
				let bond= bonds[i];
				let atoms_refs= bond.getAttribute("atomRefs2");
				let bondorder= parseFloat(bond.getAttribute("order"));

				//Find out the atoms referred in atoms_refs
				let atomsbonded_string_arr = atoms_refs.split(" ");
				if (atomsbonded_string_arr.length==2){
					//find the LPAtom objects already stored in LPMolecule.listAtoms
					
					m_molecule.addBondWithAtomsSerialNumbers(atomsbonded_string_arr[0] , atomsbonded_string_arr[1]);

				}

			}
		}
		return m_molecule;
	}

	static getLPMoleculeFromArrBuffPDBFile(arrbuff0){

		//Do not use PDBLoader because it depends on THREE.js.
		//Rewriten the code to parse the PDB file

		//Transform blob to text and then reads line-by-line

		/*let blobU8 = new Uint8Array(blob);
		let enc = new TextDecoder("utf-8"); //TODO: Edge browser does not suport TextDecoder for unknown reason.
		let dec_text = enc.decode(blobU8);
		*/

		var dec_text= this.getTextFromArrayBuffer(arrbuff0);
		
		let new_molecule=null;

		//Reads line by line
		let blob_text_Lines = dec_text.split('\n'); //Splits header to lines

		if (blob_text_Lines.length>0){
			new_molecule= new LPMolecule();

			for ( let iline = 0; iline < blob_text_Lines.length ; iline ++ ) {
				let sline = blob_text_Lines[iline];
	
				if ( sline.substr( 0, 4 ) === 'ATOM' || sline.substr( 0, 6 ) === 'HETATM' ) {
	
					let x = parseFloat( sline.substr( 30, 8 ) );
					let y = parseFloat( sline.substr( 38, 8 ) );
					let z = parseFloat( sline.substr( 46, 8 ) );
					let ser_number_string = sline.substr( 6, 5 );
	
					let pdb_elname = sline.substr( 12, 2 ).trim().toLowerCase();
	
					//Create atom here
					let at0 = new LPAtom(pdb_elname, x,y,z,ser_number_string);
					
					new_molecule.addAtom(at0);
				}
				if (sline.substr(0,6) === "CONECT" ) {
					let a1_s= sline.substr(6, 5);
					let a2_s= sline.substr(11,5);	
					new_molecule.addBondWithAtomsSerialNumbers(a1_s, a2_s);

					if (sline.length >= 21){
						let a3_s = sline.substr(16,5);
						new_molecule.addBondWithAtomsSerialNumbers(a1_s, a3_s);

						if (sline.length >= 26){
							let a4_s = sline.substr(21,5);
							new_molecule.addBondWithAtomsSerialNumbers(a1_s,a4_s);

							if(sline.length >= 31){
								let a5_s = sline.substr(26,5);
								new_molecule.addBondWithAtomsSerialNumbers(a1_s, a5_s);
							}
						}
					}
				}

				if (sline.search("COMPND") !=-1 ){
					new_molecule.moleculeName = sline.substr(7).trim();
				}
					
				//Do not process "END" or "TER" lines
			}
		}
		return new_molecule;
	}

	addBondWithAtomsSerialNumbers(sn1, sn2 , bondorder){
		//Assumes sn1 and sn2 are strings
		
		let lpat1= this.listAtoms.find(
			function( lpatom ) {
				return lpatom.ser_number_string.trim().toLowerCase() == sn1.trim().toLowerCase();
			}
		)
		let lpat2= this.listAtoms.find(
			function( lpatom ) {
				return lpatom.ser_number_string.trim().toLowerCase() == sn2.trim().toLowerCase();
			}
		)
		
		if (lpat1 != undefined && lpat2!= undefined){
			let new_bond= new LPBond(lpat1 , lpat2 , bondorder);
			this.addBond( new_bond )
		}
	}

	centerMolecule(){
		//Finds the center of the molecule and adjust atom positions to center molecule

		//Calculates average x,y,z coordinates
		var sumX=0;
		var sumY=0;
		var sumZ=0; 
		for (let i=0; i<this.listAtoms.length ; i++){
			let at0 = this.listAtoms[i];
			sumX += at0.x_Angst; 
			sumY += at0.y_Angst;
			sumZ += at0.z_Angst;


		}
		let avx = sumX / this.listAtoms.length;
		let avy = sumY / this.listAtoms.length;
		let avz = sumZ / this.listAtoms.length;

		//Uses the average xyz coordinates to correct atom positions
		for (let i=0; i<this.listAtoms.length ; i++){
			this.listAtoms[i].x_Angst -= avx ;
			this.listAtoms[i].y_Angst -= avy ;
			this.listAtoms[i].z_Angst -= avz ;
		}
	}

	getMoleculeBoundingSphereRadius(){
		//Calculates the bounding sphere with center at origin
		var rmax=0;

		for (let i=0; i<this.listAtoms.length; i++) {
			let at0 = this.listAtoms[i];

			let r0 = Math.sqrt(
				(at0.x_Angst * at0.x_Angst) + 
				(at0.y_Angst * at0.y_Angst) + 
				(at0.z_Angst * at0.z_Angst)
				)
			
				rmax= Math.max( rmax , r0);
		}

		return rmax;

	}

	deepClone(){
		//Creates a clone of LPMolecule
		//Useful when duplicating

		var cl_molecule= new LPMolecule();
		
		//Clones atoms individually
		this.listAtoms.forEach(function(at){
			cl_molecule.addAtom(at.clone());
		});

		//Clones bonds, but be careful with atom indexes
		//this.listBonds.forEach() will make 'this' undefined
		//reason why i need to use for loop
		for (let i=0; i< this.listBonds.length ; i++){
			let bond0= this.listBonds[i];
			let i1 = this.listAtoms.indexOf(bond0.atom1);
			let i2 = this.listAtoms.indexOf(bond0.atom2);
			
			//Same atoms index in new molecule
			let newbond = new LPBond(cl_molecule.listAtoms[i1] , cl_molecule.listAtoms[i2]);
			cl_molecule.addBond(newbond);
		}

		cl_molecule.moleculeName = this.moleculeName ;
		cl_molecule.filename = this.filename;

		return cl_molecule;

	}

	mirrorYZ(){
		this.listAtoms.forEach(function(at0){
			at0.x_Angst = -at0.x_Angst;
		})
	}
}
