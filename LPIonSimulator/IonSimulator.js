//classes and functions to do ion simulations
//By Luis Perdigao 2020

class CPotentialArray {

	constructor(pa_length_x = 40, pa_length_y=40 , pa_length_z = 40, dx=1.0 , dy=1.0, dz=1.0 , domElement=document){
		this.pa_length_z= pa_length_z;
		this.pa_length_y= pa_length_y;
		this.pa_length_x= pa_length_x;
		this.dx= dx;
		this.dy= dy;
		this.dz= dz;

		this.pa_max=0;
		this.pa_min=0;

		this.potential= [];
		this.isElectrode = [];

		this.needsRefining = true;

		this.convergeMaxChange=0;
		this.convergenceCriteria = 1e-3;
		this.convergeNiterations = 0;
		this.convergeMaxIterations = 100000;
		this.isConverging = false;
		this.isConverged = false;

		//not implemented, probably not needed
		//this.domElement = domElement; //dom Element to raise custom events
		//this.EPConvergeCompletedEvent= new Event("EPConvergeCompletedEvent");
		//this.EPElectrodesChangedEvent = new Event("EPElectrodesChangedEvent");

		this.initArrays();

	}

	initArrays(){
		this.potential = new Array(this.pa_length_z);
		this.isElectrode = new Array(this.pa_length_z);

		for (let iz = 0 ; iz < this.pa_length_z ; iz++ ){
			this.potential[iz] = new Array(this.pa_length_y);
			this.isElectrode[iz] = new Array(this.pa_length_y);

			for (let iy = 0 ; iy < this.pa_length_y ; iy++){
				this.potential[iz][iy] = new Array(this.pa_length_x);
				this.isElectrode[iz][iy] = new Array(this.pa_length_x);

				for (let ix=0 ; ix < this.pa_length_x ; ix++ ){	
					this.potential[iz][iy][ix] = 0;
					this.isElectrode[iz][iy][ix] = false;
				}
			}
		}
		//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that PA changed (maybe not needed here)
	}

	//Add electrode to the potential array, with limits defined by indices (ix0,iy0,iz0) - (ix1,iy1,iz1)
	//including last indexes
	//if badd= true then add electrode, otherwise, remove
	addRemoveBoxElectIndex( ix0,ix1 , iy0,iy1 , iz0,iz1 , elpotvalue, badd=true){
		
		if (ix0<0) ix0=0;
		if (ix1>=this.pa_length_x) ix1= this.pa_length_x-1;
		if (iy0<0) iy0=0;
		if (iy1>=this.pa_length_y) iy1= this.pa_length_y-1;
		if (iz0<0) iz0=0;
		if (iz1>=this.pa_length_z) iz1= this.pa_length_z-1;

		if (ix1>=ix0 && iy1>=iy0 && iz1>=iz0){
			this.needsRefining = true;

			for (let iz = iz0 ; iz<=iz1 ; iz++){
				for (let iy=iy0 ; iy<=iy1 ; iy++){
					for (let ix=ix0 ; ix<=ix1 ; ix++){
						this.potential[iz][iy][ix]= elpotvalue;
						this.isElectrode[iz][iy][ix] = badd;
					}
				}
			}
		}else{
			console.debug("CPotentialArray.addRemoveBoxElectIndex(), second index not higher than first index.")
		}

		//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that PA changed

	};

	setElectrodePotentialWithIndex(ix,iy,iz, elpotvalue, badd=true){
		if (ix>=0 &&
			ix< this.pa_length_x &&
			iy>=0 &&
			iy < this.pa_length_y &&
			iz>=0 &&
			iz< this.pa_length_z){
				this.potential[iz][iy][ix]= elpotvalue;
				this.isElectrode[iz][iy][ix] = badd;
			}
	}


	//Add box electrode to the potential array, with limits defined by positions (x0,y0,z0) - (x1,y1,z1)
	//if badd= true then add electrode, otherwise, remove
	addRemoveBoxElect( x0,x1 , y0,y1 , z0,z1 , elpotvalue, badd=true){
		//Find out corresponding indexes
		let ix0 = Math.round(x0/this.dx);
		let ix1 = Math.round(x1/this.dx);
		let iy0 = Math.round(y0/this.dy);
		let iy1 = Math.round(y1/this.dy);
		let iz0 = Math.round(z0/this.dz);
		let iz1 = Math.round(z1/this.dz);

		if (ix0<0) ix0=0;
		if (ix1>this.pa_length_x) ix1= this.pa_length_x-1;
		if (iy0<0) iy0=0;
		if (iy1>this.pa_length_y) iy1= this.pa_length_y-1;
		if (iz0<0) iz0=0;
		if (iz1>this.pa_length_z) iz1= this.pa_length_z-1;

		if (ix1>=ix0 && iy1>=iy0 && iz1>=iz0){
			this.needsRefining = true;

			for (let iz = iz0 ; iz<=iz1 ; iz++){
				for (let iy=iy0 ; iy<=iy1 ; iy++){
					for (let ix=ix0 ; ix<=ix1 ; ix++){
						this.setElectrodePotentialWithIndex(ix,iy,iz, elpotvalue, badd);
					}
				}
			}
		}else{
			console.error("CPotentialArray.addRemoveBoxElect(), second position not higher than first electrode position.");
		}
		//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that PA changed
	};

	//Add box electrode to the potential array, with limits defined by positions (x0,y0,z0) - (x1,y1,z1)
	//if badd= true then add electrode, otherwise, remove
	addRemoveCylElect( x,y,z, r,h, dir , elpotvalue, badd=true){
		if (r>0 && h>0){
			if (dir=="x" || dir=="y" || dir=="z"){
				//Find out corresponding indexes

				//let iy = Math.round(y/this.dy);
				//let iz = Math.round(z/this.dz);

				let r2 = r*r;

				this.needsRefining = true;

				if (dir=="x"){
					//x direction will be the height direction
					let ix0 = Math.round(x/this.dx);
					let ix1 = Math.round((x+h)/this.dx);
					//TODO for loop on x indexes.
					
					let y0 = y-r; //absolute poistion of start and end of circle
					let y1 = y+r;
					
					let iy0 = Math.round(y0/this.dy);//abs position in index units
					let iy1 = Math.round(y1/this.dy);

					for (let jy= iy0; jy <= iy1 ; jy++){
						let dy2= (jy*this.dy-y) ** 2; //relat to centre position squared
						//get limits on the other axis
						let dz1 =  Math.sqrt( r2 - dy2);

						let iz0 = Math.round((z-dz1)/this.dz);
						let iz1 = Math.round((z+dz1)/this.dz);
						
						for (let jz = iz0 ; jz <= iz1 ; jz++){
							for (let jx = ix0 ; jx <= ix1 ; jx++){
								this.setElectrodePotentialWithIndex(jx,jy,jz, elpotvalue, badd);				
							}
						}
					}
				}else if( dir=="y" ){
					//y direction will be the height direction
					let iy0 = Math.round(y/this.dy);
					let iy1 = Math.round((y+h)/this.dy);
					
					let x0 = x-r; //absolute position of start and end of circle
					let x1 = x+r;
					
					let ix0 = Math.round(x0/this.dx);//abs position in index units
					let ix1 = Math.round(x1/this.dx);

					for (let jx= ix0; jx <= ix1 ; jx++){
						let dx2= (jx*this.dx-x) ** 2; //relat to centre position squared
						//get limits on the other axis
						let dz1 =  Math.sqrt( r2 - dx2);

						let iz0 = Math.round((z-dz1)/this.dz);
						let iz1 = Math.round((z+dz1)/this.dz);
						
						for (let jz = iz0 ; jz <= iz1 ; jz++){
							for (let jy = iy0 ; jy <= iy1 ; jy++){
								this.setElectrodePotentialWithIndex(jx,jy,jz, elpotvalue, badd);				
							}
						}
					}
				}else if( dir=="z" ){
					//y direction will be the height direction
					let iz0 = Math.round(z/this.dz);
					let iz1 = Math.round((z+h)/this.dz);
					
					let x0 = x-r; //absolute position of start and end of circle
					let x1 = x+r;
					
					let ix0 = Math.round(x0/this.dx);//abs position in index units
					let ix1 = Math.round(x1/this.dx);

					for (let jx= ix0; jx <= ix1 ; jx++){
						let dx2= (jx*this.dx-x) ** 2; //relat to centre position squared
						//get limits on the other axis
						let dy1 =  Math.sqrt( r2 - dx2);

						let iy0 = Math.round((y-dy1)/this.dy);
						let iy1 = Math.round((y+dy1)/this.dy);
						
						for (let jy = iy0 ; jy <= iy1 ; jy++){
							for (let jz = iz0 ; jz <= iz1 ; jz++){
								this.setElectrodePotentialWithIndex(jx,jy,jz, elpotvalue, badd);
							}
						}
					}
				}
				//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that PA changed

			}else{
				console.error("Cylinder direction not understood.")
			}
		}else{
			console.error("Radius r or height h of cylinder are not positive values.");
		}

	};

	convergeIteration() {
		//added support for arrays with different dx dy and dz
		
		let maxchange = 0;

		//TODO: prevent division by zero. In javascript will lead to inf values.
		let oneoverdx2 = 1/this.dx/this.dx ;
		let oneoverdy2 = 1/this.dy/this.dy ;
		let oneoverdz2 = 1/this.dz/this.dz ;

		for (let iz = 0 ; iz < this.pa_length_z ; iz++){
			for (let iy = 0 ; iy < this.pa_length_y ; iy++){
				for (let ix = 0 ; ix < this.pa_length_x ; ix++){
					let numerator=0;
					let denominator=0;
					let newval=0;
					
					if ( ! this.isElectrode[iz][iy][ix]  ){
						if (iz>0){
							numerator += this.potential[iz-1][iy][ix] * oneoverdz2 ;
							denominator += oneoverdz2;
						}
						if (iz< this.pa_length_z-1){
							numerator += this.potential[iz+1][iy][ix] * oneoverdz2;
							denominator += oneoverdz2 ;						
						}
						if (iy>0){
							numerator += this.potential[iz][iy-1][ix] * oneoverdy2;
							denominator += oneoverdy2 ;
						}
						if (iy< this.pa_length_y-1){
							numerator += this.potential[iz][iy+1][ix] * oneoverdy2;
							denominator += oneoverdy2 ;						
						} 
						if (ix>0){
							numerator += this.potential[iz][iy][ix-1] * oneoverdx2;
							denominator += oneoverdx2;
						}
						if (ix< this.pa_length_x-1){
							numerator += this.potential[iz][iy][ix+1] * oneoverdx2;
							denominator += oneoverdx2;						
						}

						if (denominator>0){
							newval = numerator / denominator ;
							let vchange = Math.abs( newval - this.potential[iz][iy][ix] ) ;
							maxchange = Math.max( vchange , maxchange);
							this.potential[iz][iy][ix] = newval;
						}
					}
				}
			}
		}
		this.convergeMaxChange = maxchange;
	};


	converge(){
		//Consider making this function a 'promise'

		//while loop written as a for loop
		this.convergeIteration();
		this.convergeNiterations++;
		this.isConverging=true;
		this.isConverged=false; //default

		if (this.convergeMaxChange <= this.convergenceCriteria){
			this.isConverging=false; //Signals that calculation of EP is complete
			this.isConverged=true;
			console.log("converge() Convergence criteria reached.");
			this.updateMaxMin();
			//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that EP changed
			this.domElement.dispatchEvent(this.EPConvergeCompletedEvent);
		}else if (this.convergeNiterations > this.convergeMaxIterations ){
			console.log("converge() Reached maximum iterations. Stopping refining.");
			this.isConverging=false;
			this.updateMaxMin();
			//this.domElement.dispatchEvent(this.EPElectrodesChangedEvent); //signals that EP changed
		}else{
			//setTimeout ( this.converge , 10 ) ; //10ms between each iteration to give some time for graphics to update
			//setTimeout causes javascript to use the global scope, so it will fail.
			setTimeout ( this.converge.bind(this) , 0 ) ; //0ms between each iteration to give some time for graphics to update
		}
	}

	updateMaxMin(){
		//First point resets values
		let v0 = this.potential[0][0][0];

		this.pa_min=v0;
		this.pa_max=v0

		for (let iz = 0 ; iz < this.pa_length_z ; iz++){
			for (let iy = 0 ; iy < this.pa_length_y ; iy++){
				for (let ix = 0 ; ix < this.pa_length_x ; ix++){
					v0=this.potential[iz][iy][ix];
					this.pa_max = Math.max(v0, this.pa_max);
					this.pa_min = Math.min(v0, this.pa_min);
				}
			}
		}
	}

	isParticlePositionValid(x,y,z){
		let isValid = true; //default is valid, unless...

		//Check if point is inside limits of potential array
		if ( x<0 || x> this.pa_length_x * this.dx*this){
			isValid = false;
		}
		if ( y<0 || y> this.pa_length_y * this.dy*this){
			isValid = false;
		}
		if ( z<0 || z> this.pa_length_z * this.dz*this){
			isValid = false;
		}

		//check if point is inside electrode
		if ( this.isElectrode[ Math.floor(z / this.dz) ] [ Math.floor( y / this.dy) ][ Math.floor(this.dx) ]){
			isValid=false;
		}

		return isValid;
	}

	getEfieldAt(x,y,z){
		let Ex=0;
		let Ey=0;
		let Ez=0;

		let isOutsidePA=false;
		let isInsideElectrode = false;

		//Convert position to 'index' coordinates
		let ixfloat = x / this.dx;
		let iyfloat = y / this.dy;
		let izfloat = z / this.dz;

		if (ixfloat < -0.5 || iyfloat < -0.5 || izfloat < -0.5 ||
			ixfloat >= (this.pa_length_x-0.5) || iyfloat>= (this.pa_length_y-0.5) || izfloat >= (this.pa_length_z-0.5) ) {
				isOutsidePA = true;
				//console.debug("point is outside PA");
				// ... and exit function
		}else{
			//Check it is inside electrode
			
			//origin index, closest index point
			let ix0 = Math.round( ixfloat );
			let iy0 = Math.round( iyfloat );
			let iz0 = Math.round( izfloat );
			
			if (this.isElectrode[iz0][iy0][ix0]){
				isInsideElectrode=true;
				//console.debug("point is inside electrode");
			}else{
				//OK to calculate Ex Ey and Ez
				if (ix0==0 || ix0 == this.pa_length_x-1){
					Ex=0;
				}else{
					let ix1=ix0-1; 
					//Check for a point in the same direction
					//TODO: This is a very crude method
					if (ixfloat>=ix0){
						ix1=ix0+1;
					}
					Ex = (this.potential[iz0][iy0][ix0] - this.potential[iz0][iy0][ix1]) / (this.dx * (ix1-ix0)) ;
				}
				if (iy0==0 || iy0 == this.pa_length_y-1){
					Ey=0;
				}else{
					let iy1=iy0-1; 
					//Check for a point in the same direction
					//TODO: This is a very crude method, but probably ok
					if (iyfloat >= iy0){
						iy1=iy0+1;
					}
					Ey = (this.potential[iz0][iy0][ix0] - this.potential[iz0][iy1][ix0]) / (this.dy * (iy1-iy0)) ;
				}
				if (iz0==0 || iz0 == this.pa_length_z-1){
					Ez=0;
				}else{
					let iz1=iz0-1; 
					//Check for a point in the same direction
					//TODO: This is a very crude method
					if (izfloat >= iz0){
						iz1=iz0+1;
					}
					Ez = (this.potential[iz0][iy0][ix0] - this.potential[iz1][iy0][ix0]) / (this.dz * (iz1-iz0)) ;
				}
			}
		}

		return {Ex , Ey ,  Ez, isInsideElectrode , isOutsidePA};
	}

	static newFromCPotentialArray(cparray){
		//Creates a new CPotentialArray object from another CPotentialArray
		//This is a shallow copy but it is very useful when deserializing objects
		//since when serialising/deserialising, the methods are lost

		var newcparray = new CPotentialArray(0,0,0)//Initialise empty

		//Copy all properties (based in constructor)
		newcparray.pa_length_z= cparray.pa_length_z;
		newcparray.pa_length_y= cparray.pa_length_y;
		newcparray.pa_length_x= cparray.pa_length_x;
		newcparray.dx= cparray.dx;
		newcparray.dy= cparray.dy;
		newcparray.dz= cparray.dz;


		newcparray.potential= cparray.potential;
		newcparray.isElectrode = cparray.isElectrode;

		newcparray.pa_max= cparray.pa_max;
		newcparray.pa_min= cparray.pa_min;

		newcparray.needsRefining = cparray.needsRefining;

		newcparray.convergeMaxChange=cparray.convergeMaxChange;
		newcparray.convergenceCriteria = cparray.convergenceCriteriay;
		newcparray.convergeNiterations = cparray.convergeNiterations;
		newcparray.convergeMaxIterations = cparray.convergeMaxIterations;
		newcparray.isConverging = cparray.isConverging;
		newcparray.isConverged = cparray.isConverged;


		//not implemented, probably not needed
		//newcparray.domElement = cparray.domElement; //dom Element to raise custom events
		//newcparray.EPConvergeCompletedEvent= new Event("EPConvergeCompletedEvent");
		//newcparray.EPElectrodesChangedEvent = new Event("EPElectrodesChangedEvent");

		return newcparray;
	}

} //class CPotentialArray

class CIonSimulator{
	constructor(pa_array_obj=null){
		//Initialise the Class, needs a CPotentialArray 
		this.EP=pa_array_obj;
		
		this.particles = [];
		this.timeStep_s = 1e-2 ; //10 ms

		this.timeElapsed_s = 0;
		this.iteration = 0;
		this.timeLastIteration=null;

		this.particleSources=[]; //Particles starting positions

		this.flyCount=0; //in case of continous fly, this counts the number of flys so far

		this.bSimulationStop=false; //used to flag to stop continous simulation
		this.bSimulationAutoRestart=true;
		this.bSimulationIsRunning = false;

		this.turbo=true;

		//this.integrationMethod=""; //TODO: Maybe is faster using callback function.
	}

	addParticleSource( CParticle_obj ){
		//Check type of variable
		if (CParticle_obj.constructor.name =="CParticle"){
			this.particleSources.push(CParticle_obj);

			let pcopy = CParticle_obj.clone();
			//Object.assign( {} , CParticle_obj); //Shallow copy not good enough.

			this.particles.push(pcopy);
		}
	}

	clearParticleSources(){
		this.particleSources = [];
		this.particles = [];
	}

	restartParticlesFromSources(){
		for (let i = 0 ; i < this.particleSources.length ; i++){
			this.particles[i].copy(this.particleSources[i]);
		}
		this.timeElapsed_s=0;
		this.iteration=0;
	}

	
	//Runs a single time step iteration on all particles
	runiterationstep(){

		//Uses this.timeStep_s to determine next positions

		//Each particle,
		this.particles.forEach( (pobj) => {

			//Get position.
			let x0= pobj.x_m;
			let y0= pobj.y_m;
			let z0= pobj.z_m;


			//Get the force F = q . - grad V = q.E
			let efield = this.EP.getEfieldAt(x0,y0,z0);

			if (!efield.isInsideElectrode && !efield.isOutsidePA){
				//Calculate new positions

				//Simple Euler-Crommer method
				//Apply force, adjust velocity
				let temp0 = pobj.charge_C / pobj.mass_kg * this.timeStep_s;

				//x
				let v1x = pobj.vx_ms + temp0 * efield.Ex ;
				let x1 =  pobj.x_m + this.timeStep_s * v1x ;

				//y
				let v1y = pobj.vy_ms + temp0 * efield.Ey ;
				let y1 =  pobj.y_m + this.timeStep_s * v1y ;

				//z
				let v1z = pobj.vz_ms + temp0 * efield.Ez ;
				let z1 =  pobj.z_m + this.timeStep_s * v1z ;

				//Set the new values
				pobj.x_m = x1;
				pobj.vx_ms = v1x;

				pobj.y_m = y1;
				pobj.vy_ms = v1y;

				pobj.z_m = z1;
				pobj.vz_ms = v1z;

			}else{
				//particle is inside electrode or outsidePA
				if (efield.isInsideElectrode){
					pobj.status = 1;
				}else{
					pobj.status = 2;
				}
			}
			//next Particle
		});
		this.iteration++;
		this.timeElapsed_s += this.timeStep_s;
	}

	//Run iterations continuously until all particles hit boundaries or electrodes
	// bAutoStart is to specify whether to restart all particles from sources after all already hit
	startSimulation(){
		if (this.bSimulationIsRunning == false){ //Don't start a simulation while one is running
			this.bSimulationStop=false;
			this.iteration=0;
			this.timeElapsed_s=0;
			this.cycleSimulation();
		}


	}
	cycleSimulation(){
		this.bSimulationIsRunning=true;
		if (!this.bSimulationStop ){
			let t0 = performance.now();
			let tms = 0;
			let t1 = 0 ;

			do{
				this.runiterationstep();
				t1= performance.now();
				tms = t1-t0;
			} while (this.turbo ==true && tms<16)
			//every 16ms let the javascript main loop run for redraw graphics

			//Tries to find out if all the particles hit
			let bAllParticlesHit= true;
			this.particles.forEach( p0 => {
				if (p0.status == CParticle.statusvalues.OK){
					bAllParticlesHit=false;
				}
			});
	
			//If all particles hit and bAutostart enabled then restart all particles
			if (bAllParticlesHit && this.bSimulationAutoRestart){
				this.restartParticlesFromSources();
			}

			setTimeout( this.cycleSimulation.bind(this) , 0 ); //runs continously, but does not hang main loop
			//Note the use of bind(this). It is important when using setTimeout inside a class
		}else{
			console.debug( "cycleSimulation(): Stop signal recieved.")
			this.bSimulationIsRunning=false;
		}
	}
	stopSimulation(){
		this.bSimulationStop=true;
	}

	static newFromCIonSimulator(cionsim){
		//Creates a new CIonSimulator from the one given
		//Used when deserialising, to restore methods

		var newionsim = new CIonSimulator(); 

		//Copy the Potential Array (Electric potential EP)
		newionsim.EP= CPotentialArray.newFromCPotentialArray(cionsim.EP);

		//particles, add them one by one
		newionsim.particleSources=[];
		cionsim.particleSources.forEach( function(ps0){
			var newparticle = CParticle.newFromCParticle(ps0);
			newionsim.addParticleSource(newparticle);
		} )

		newionsim.timeStep_s = cionsim.timeStep_s ; //10 ms
		newionsim.timeElapsed_s = 0;
		newionsim.iteration = 0;
		newionsim.timeLastIteration=null;
		newionsim.flyCount=0; //in case of continous fly, newionsim counts the number of flys so far
		newionsim.bSimulationStop=false; //used to flag to stop continous simulation
		newionsim.bSimulationAutoRestart=true;
		newionsim.bSimulationIsRunning = false;
		newionsim.turbo=cionsim.turbo;

		return newionsim;
	}
}

//Particle as a class
class CParticle{
	static statusvalues={
		OK: 0,
		OUTSIDEPA: 1,
		INSIDEELECTRODE: 2
	};
	constructor(x_m ,y_m , z_m , vx_ms,vy_ms, vz_ms, charge_C , mass_kg){
		this.x_m = x_m;
		this.y_m = y_m;
		this.z_m = z_m;
		this.vx_ms = vx_ms;
		this.vy_ms = vy_ms;
		this.vz_ms = vz_ms;
		this.charge_C = charge_C;
		this.mass_kg = mass_kg;

		this.status = CParticle.statusvalues.OK;
	};

	KineticEnergy(){
		return ( (this.vx_ms*this.vx_ms +  this.vy_ms*this.vy_ms + this.vz_ms*this.vz_ms) * this.mass_kg * 0.5 );
	}

	clone(){
		return new CParticle(
			this.x_m,
			this.y_m ,
			this.z_m ,
			this.vx_ms ,
			this.vy_ms ,
			this.vz_ms ,
			this.charge_C ,
			this.mass_kg
		)
	}
	
	//copies properties of CParticle parameter to this object
	copy(cPartSrc){
		this.x_m = cPartSrc.x_m;
		this.y_m = cPartSrc.y_m;
		this.z_m = cPartSrc.z_m;
		this.vx_ms = cPartSrc.vx_ms;
		this.vy_ms = cPartSrc.vy_ms;
		this.vz_ms = cPartSrc.vz_ms;
		this.charge_C = cPartSrc.charge_C;
		this.mass_kg = cPartSrc.mass_kg;
		this.status = cPartSrc.status;
	}

	static newFromCParticle(cpart){
		var newpart = new CParticle( cpart.x_m, cpart.y_m, cpart.z_m,
			cpart.vx_ms, cpart.vy_ms, cpart.vz_ms, cpart.charge_C, cpart.mass_kg);
		newpart.status = cpart.status;

		return newpart;
	}
}
