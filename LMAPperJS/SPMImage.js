class SPMImage{
	constructor(){
		this.xsize_nm=0;
		this.ysize_nm=0;
		this.xpixels=0;
		this.ypixels=0;
		this.name="";
		this.spmImageZData=[];
		//this.imageBMP=[];
		this.imagedata_max=0;
		this.imagedata_min=0;
		this.filename="";
		this.isBitmap = false;
		this.imageAsCanvas = []; //Try storing image as canvas rather than imageBMP
	}
	
	//SPM image processing routines
	
	globalPlane() {
        //Correct image using Global plane
        //Calculates fit plane using Crammer's rule
        let SX2 = 0,
          SXY = 0, SXZ = 0,
          SY2 = 0,
          SYZ = 0,
          SZ = 0;

        let SX=0, SY=0;

        for (let iy = 0; iy < this.ypixels; iy++) {
          for (let ix = 0; ix < this.xpixels; ix++) {
            let z = this.spmImageZData[iy][ix];
            SX2 += ix * ix;
            SXY += ix * iy;
            SXZ += ix * z;
            SY2 += iy * iy;
            SYZ += iy * z;
            SZ += z;
            SX += ix;
            SY += iy;
          }
        }

        //Use Crammer to get plane fitting parameters
        //Create matrices
        const X = [[SX2, SXY, SX],
        [SXY, SY2, SY],
        [SX, SY, this.xpixels * this.ypixels]
        ];

        
        const Xt=math.transpose(X);

        //var detX = math.det(X);
        const detX= determinant3x3(X);

        if (detX == 0) {
          console.log("detX=0, cannot continue doing Crammer rule to determine best fit plane");
        } else {
          //Use Crammer with SXZ, SYZ, SZ
          var Y = [SXZ, SYZ, SZ];

          /*
          console.log("System of equation to solve is");
          for (let i0=0; i0<3; i0++){
            let Xl=X[i0];
            console.log(Xl[0].toString() + " x + "+ Xl[1].toString() + " y + " + Xl[2].toString() + " z = "+ Y[i0]);
          }*/
          
          //Makes copies of arrays and modifies relevant column modification
          var A0 = Xt.slice();
          A0[0] = Y;
          var A1 = Xt.slice();
          A1[1] = Y;
          var A2 = Xt.slice();
          A2[2] = Y;
          
          /*
          var detA0 = math.det(A0);
          var detA1 = math.det(A1);
          var detA2 = math.det(A2);
          */

          const detA0=determinant3x3(A0);
          const detA1=determinant3x3(A1);
          const detA2=determinant3x3(A2);

          //Define object Plane p = ax + by + c
          const plane_a= detA0 / detX, plane_b = detA1 / detX , plane_c= detA2 / detX;

          //var plane0 = { a: detA0 / detX, b: detA1 / detX, c: detA2 / detX };

          //console.log("Best fitted plane y = ax + by + c; plane0= " + plane0.toString());
          console.log("Best fitted plane y = ax + by + c; a=" + plane_a.toString() + " b="+plane_b.toString() + " c="+plane_c.toString() );

          //Subtract plane, and simultaneously gets max and min
          //var imagedata0 = this.imagedata.slice(); //copy image
          var imagedata0min = 0; imagedata0max = 0;
          for (let iy = 0; iy < this.ypixels; iy++) {
            for (let ix = 0; ix < this.xpixels; ix++) {
              let vplane= ((iy * plane_b) + (ix * plane_a) + plane_c);
              let value = this.spmImageZData[iy][ix] - vplane;
              //imagedata0[iy][ix] = value;
              this.spmImageZData[iy][ix] = value; //new value
              //this.imagedata[iy][ix]-=plane_c; //try

              if (iy == 0 && ix == 0) {
                imagedata0max = value;
                imagedata0min = value;
              } else {
                imagedata0max = Math.max(imagedata0max, value);
                imagedata0min = Math.min(imagedata0min, value)
              }
            }
          } //By here, image0 should have the 'planned' data, and the max and min values determined
          console.log("imagedata0max = " + imagedata0max.toString() + " ; imagedata0min = " + imagedata0min.toString());

          //Update this object
          //this.imagedata = imagedata0;
          this.imagedata_max = imagedata0max;
          this.imagedata_min = imagedata0min;

          this.updateBMP(); //Updates BMP
        }
    }
	
	flatten(){
        //Corrects line-by-line
        var SX = this.xpixels * (this.xpixels - 1) / 2; //Pythagoras formula (single line)

        var imagedata0min = 0, imagedata0max = 0;

        for (let iy = 0; iy < this.ypixels; iy++) {
          var SZ = 0,
            SX2 = 0,
            SZ2 = 0,
            SXZ = 0;
          for (let ix = 0; ix < this.xpixels; ix++) {
            var z = this.spmImageZData[iy][ix];
            SZ += z;
            SX2 += ix * ix;
            SZ2 += z * z;
            SXZ += ix * z;
          }

          /*
          Dim line0 As Line
            Dim d0 As Double = w * SX2.sum - SX.sum ^ 2
            line0.c = (SY.sum * SX2.sum - SX.sum * SXY.sum) / d0
            line0.m = (w * SXY.sum - SX.sum * SY.sum) / d0
            */
          let d0 = this.xpixels * SX2 - SX*SX;
          let line_c= (SZ * SX2 - SX * SXZ) / d0;
          let line_m= (this.xpixels * SXZ - SX * SZ) / d0;

          //Subtract line here, and calculates limits
          for (let ix = 0; ix < this.xpixels; ix++) {
            let z = this.spmImageZData[iy][ix];
            let newvalue = z - line_m * ix - line_c;
            this.spmImageZData[iy][ix] = newvalue;

            if (iy == 0 && ix == 0) {
              imagedata0max = newvalue;
              imagedata0min = newvalue;
            } else {
              imagedata0max = Math.max(imagedata0max, newvalue);
              imagedata0min = Math.min(imagedata0min, newvalue);
            }
          }
          this.imagedata_max = imagedata0max;
          this.imagedata_min = imagedata0min;
        }
        this.updateBMP(); //Updates BMP
    }

	updateBMP(){
		if (this.xpixels>0 && this.ypixels>0 && this.spmImageZData!=null){
        //   this.imageBMP = new ImageData(this.xpixels, this.ypixels);
        //   var imageRGB = this.imageBMP; //Gets the reference

		var imageRGB = new ImageData(this.xpixels, this.ypixels);

          var deltamaxmin = this.imagedata_max - this.imagedata_min;
          if (deltamaxmin == 0) {
            console.log("Error, deltamaxmin = 0, cannot create a suitable grayscale image")
          } else {
            for (let iy = 0; iy < this.ypixels; iy++) {
              for (let ix = 0; ix < this.xpixels; ix++) {
                let value = this.spmImageZData[iy][ix];
                let value256 = Math.floor((value - this.imagedata_min) / deltamaxmin * 256);

                //direct pixel manipulation
                let index = (ix + iy * this.xpixels) * 4;
                imageRGB.data[index] = value256;
                imageRGB.data[index + 1] = value256;
                imageRGB.data[index + 2] = value256;
                imageRGB.data[index + 3] = 255;
              }
            }
			console.log("updateBMP, imageRGB.data created");
			
			//Create (invisible) canvas element this.imageAsCanvas
			let newCanvas = document.createElement('canvas');
			newCanvas.width = imageRGB.width;
			newCanvas.height = imageRGB.height;
	
			let ctx = newCanvas.getContext("2d");
			ctx.putImageData(imageRGB, 0, 0);
			
			this.imageAsCanvas=newCanvas;
			
			this.isBitmap=false;
	

          }
        }else{
          console.debug("Cannot update BMP");
        }
    }

    determinant3x3(arr){
      //Assumes that the input variable is a 3x3 array
      return arr[0][0] * (arr[1][1]*arr[2][2] - arr[1][2]*arr[2][1]) -
                  arr[0][1] * (arr[1][0]*arr[2][2] - arr[1][2]*arr[2][0]) +
                  arr[0][2] * (arr[1][0]*arr[2][1] - arr[1][1]*arr[2][0]);
    }
	
	static getSPMImageFromArrBuffSXMFile(arrbuff0){
		
		var m_spmimage=new SPMImage();
		
		//arrbuff0 is a ArrayBuffer object, and cannot be accessed directly.
		//To access it, we need to say what type it is

		var arrbuff0U8 = new Uint8Array(arrbuff0);
		//try tis method to locate Chr(&H1A) & Chr(&H4), which indicates begining of binary data
		var binpos0 = arrbuff0U8.findIndex(
			function (element, index, array){
			  if ( element==26 && array[index+1]==4 ){
				return true;
			  }
			  else{
				return false;
			  }
			}
	    );
    
		if (binpos0>0){
			//Gets the header
			//var sHeaderU8 = arrbuff0U8.slice(0, binpos0);
			var sHeaderU8 = arrbuff0.slice(0, binpos0);

			/*
			var enc = new TextDecoder("utf-8");
			var dec_sHeaderU8 = enc.decode(sHeaderU8);
			*/
			var dec_sHeaderU8= this.getTextFromArrayBuffer(sHeaderU8)

			//Reads line by line
			var spHeaderLines = dec_sHeaderU8.split('\n'); //Splits header to lines
			var s="";
			//Each value starts with a title in format :<title>:
			//Followed by value or values
			var s1="";
			var s1_split="";
			var nchannels=0;
			for(var i=0; i<spHeaderLines.length ; i++){
				s=spHeaderLines[i];
				if (s.length > 1) {
					if (s.includes("SCAN_PIXELS")) {
						//next line has the number of pixels X and Y
						i+=1;
						s1 = spHeaderLines[i];
						s1_split= s1.split(/[ ]+/); //use regex, one or more spaces

						m_spmimage.xpixels=Number(s1_split[1]); //Zero index is empty, number appears at index 1
						m_spmimage.ypixels=Number(s1_split[2]);
					}else{
						if (s.includes("SCAN_RANGE")) {
							//next line has the size in X and Y
							i+=1;
							s1 = spHeaderLines[i];
							s1_split= s1.split(/[ ]+/); //use regex, one or more spaces
							m_spmimage.xsize_nm=Number(s1_split[1]) * 1e9;
							m_spmimage.ysize_nm=Number(s1_split[2]) * 1e9;
						}else{
							if (s.includes("DATA_INFO")) {
								//First line after this contains a header.
								//Check column position for Name and Direction
								//If direction says both, then it has two images.
								i+=1;
								s1 = spHeaderLines[i];
								s1_split= s1.split(/[ ]+/); //use regex, one or more spaces
								var tabDirIndex = s1_split.indexOf('Direction');
								var nameIndex=s1_split.indexOf('Direction');

								//Read following lines. If line is empty, finish the channel count
								for ( ; s1 != "" && i<spHeaderLines.length;) {
									i++;
									s1 = spHeaderLines[i];
									if (s1!=""){
										nchannels+=1;
										//s1_split= s1.split(/[ ]+/); //use regex, one or more spaces
										if (s1.indexOf('both')>=0){
											nchannels+=1;
										}
									}
								}
							}
						}
					}
				}
			}

			console.log("xpixels =" + m_spmimage.xpixels.toString() +
				", ypixels = "+m_spmimage.ypixels.toString() +
				", xsize = "+ m_spmimage.xsize_nm.toString() +
				", ysize = "+ m_spmimage.ysize_nm.toString());
			console.log("nchannels=" + nchannels.toString());

			//Check header data is ok before preceeding
			if (m_spmimage.xpixels>0 && m_spmimage.ypixels>0 &&
				m_spmimage.xsize_nm>0 && m_spmimage.ysize_nm>0 &&
				nchannels>0){
				//Read binary data
				//Need to select one of the channels, if applicable
				binpos0+=2; //reposition pointer

				//Read first channel only
				var ichannel=0;
				var channel_bytelength = m_spmimage.xpixels * m_spmimage.ypixels *4 ;
				//var data0_U8 = arrbuff0U8.slice(binpos0, channel_bytelength);

				//Convert to float32
				//const data0_F32 = Float32Array.from(data0_U8);
				const data0_DV= new DataView(arrbuff0); //try

				//Fills m_spmimage.imagedata[][] with data
				var image=new Array(m_spmimage.ypixels);
				for (let iy = 0; iy < m_spmimage.ypixels; iy++) {
					image[iy] = new Array(m_spmimage.xpixels);
					for (let ix = 0; ix < m_spmimage.xpixels; ix++) {
						//Data starts at position 3
						//Note also that iy Flips vertically by using m_spmimage.ypixels-iy-1
						let pos0 = 4*(ichannel * (m_spmimage.xpixels * m_spmimage.ypixels) + (m_spmimage.ypixels-iy-1) * m_spmimage.xpixels + ix) + binpos0;

						//var value = data0_F32[pos0];
						let value=data0_DV.getFloat32(pos0,false);
						//image[iy][ix] = value; //first element of single element array after conversion
						image[iy][ix] = value; 
					}
				}
			
				console.log("image[iy][ix] created");
				
				m_spmimage.spmImageZData=image;

				m_spmimage.flatten();
				console.log("m_spmimage flattened");
		
			}else{
				console.log("Error with parameters read in header. Setting m_spmimage.imagedata=null");
				m_spmimage=null;
			}

		}else{
			console.log("getSPMImageFromArrBuffSXMFile: Error in arrbuff0, cannot find characters &H1A &H04. No binary data. Setting m_spmimage.imagedata=null");
			m_spmimage=null;
		}

		return m_spmimage;
	}

	static getSPMImageFromArrBuffCreatecFile(arrbuff0){
	
			//TODO
			//Createc files have compressed binary in zip, so we need to use Pako library to decompress
	}

	static getTextFromArrayBuffer(arrbuff){
		let s= String.fromCharCode.apply( null, new Uint8Array(arrbuff)); // .apply() is a function that works with arrays
		return s;
	}

	static getSPMImageFromArrBuffWSXMSTPFile(arrbuff0){
		
		var m_spmimage=null;

		//WSxM stp files consist of an header, with several parameters
		//followed by binary data.

		//Binary data comes after the string "[Header end]CRLF"
		//Another way to spot the location of binary data is to read line that says
		//"Image header size: <number of bytes>"

		let headerShortArrBuff= arrbuff0.slice(0, 200); //200 should be enough
		let shortHeader= this.getTextFromArrayBuffer(headerShortArrBuff);

		//Try to locate line that says "Image header size:"
		// using regular expressions
		let re = /Image header size:\s*[0-9]+/ ; //Regular expression
		//Match string "Image header size:", followed by zero or more white space characters, followed by one or more digits.
		let imheaderString = shortHeader.match(re);

		//If match was found then extract the value provided
		re = /[0-9]+/ ;
		let imheaderHeaderSize = parseInt(shortHeader.match(re));

		if (imheaderHeaderSize>100){ // Attention! Arbitrary header size minimum
			
			//Seperate binary data from header
			let sHeaderAB =  arrbuff0.slice(0, imheaderHeaderSize);
			//Convert header to string
			let sHeader = this.getTextFromArrayBuffer(sHeaderAB);

			//Processes information, line by line
			let spHeaderLines = sHeader.split('\n'); //Splits header to lines

			var xampl=0; //in nm
			var yampl=0;
			var ncols=0;
			var nrows=0;
			var zampl=0;

			spHeaderLines.forEach( function(s){
				if (s.includes("X Amplitude:") ){
					//Gets the value after colon using regex
					let svalue =s.split(/\s*:\s*/)[1].match(/\d*\.?\d*/); //regex,colon, and then decimal
					xampl = parseFloat( svalue ); //Assumes values are in nm
					if (s.includes("Å")) xampl/=10;
				}else{
					if (s.includes("Y Amplitude:") ){
						let svalue =s.split(/\s*:\s*/)[1].match(/\d*\.?\d*/); //regex,colon, and then decimal
						yampl = parseFloat( svalue ); //Assumes values are in nm
						if (s.includes("Å")) yampl/=10;
					}else{
						if (s.includes("Number of columns:") ){
							let svalue =s.split(/\s*:\s*/)[1].match(/\d*\.?\d*/); //regex,colon, and then decimal
							ncols = parseInt( svalue ) ;
						}else{
							if (s.includes("Number of rows:") ){
								let svalue =s.split(/\s*:\s*/)[1].match(/\d*\.?\d*/); //regex,colon, and then decimal
								nrows = parseInt( svalue ) ;
							}else{
								if (s.includes("Z Amplitude:") ){
									let svalue =s.split(/\s*:\s*/)[1].match(/\d*\.?\d*/); //regex,colon, and then decimal
									zampl = parseFloat( svalue );
									if (s.includes("Å")) zampl/=10;
								}
							}
						}
					}
				}
			});

			//Check is ok to proceed
			if (ncols>0 && nrows>0 && xampl>0 && yampl>0){
				m_spmimage=new SPMImage();

				m_spmimage.xpixels=ncols;
				m_spmimage.ypixels=nrows;
				m_spmimage.xsize_nm =xampl;
				m_spmimage.ysize_nm =yampl;

				let binLengthByteExpected = ncols*nrows*4; //4bytes, 32 bit, double precision
				
				//Use DataView to extract data, but need byte position
				let binpos0=imheaderHeaderSize;
				
				const data0_DV= new DataView(arrbuff0);

				//Fills m_spmimage.imagedata[][] with data
				var image=new Array(m_spmimage.ypixels);
				for (let iy = 0; iy < m_spmimage.ypixels; iy++) {
					image[iy] = new Array(m_spmimage.xpixels);
					for (let ix = 0; ix < m_spmimage.xpixels; ix++) {
						
						//let pos0 = 8*(iy * m_spmimage.xpixels + ix) + binpos0;
						let pos0 = 8*( (m_spmimage.ypixels-iy-1) * m_spmimage.xpixels + (m_spmimage.xpixels-ix-1) )+ binpos0;

						//Lmapper SPMImageWSxMSTP.vb suggests it's a 64bit float (double precision, 8 bytes) little endian
						let value=data0_DV.getFloat64(pos0,true); 
						//image[iy][ix] = value; //first element of single element array after conversion
						image[iy][ix] = value; 
					}
				}
				m_spmimage.spmImageZData=image;
				m_spmimage.flatten();

			}else{
				console.log("Error with parameters read in header.");
			}


		}else{
			console.log("Error: Header too small");
		}

		return m_spmimage;

	}

	//Not working
	static async getSPMImageFromBitmapFile(base64EncString){
		// let m_spmimage = new SPMImage();

		// let newImage = new Image();
		// let newCanvas=document.createElement('canvas');

		// const awaitfunction = new Promise( resolve  => {
		// 		console.log("new Promise");
		// 		newImage.onload = function (){
		// 			console.log("newImage.onload");
		// 			newCanvas.width = newImage.width;
		// 			newCanvas.height = newImage.height;
		// 			let ctx = newCanvas.getContext('2d');
		// 			ctx.drawImage(newImage, 0, 0);

		// 			//No need to fill .imageBMP with data
		// 			m_spmimage.xpixels = newImage.width;
		// 			m_spmimage.ypixels = newImage.height;
		// 			m_spmimage.xsize_nm=10 ;
		// 			m_spmimage.ysize_nm=10;
		// 			m_spmimage.isBitmap=true;
			
		// 			m_spmimage.imageAsCanvas = newCanvas;

		// 			resolve();
		// 		};

		// 		console.log("pre newImage.src");
		// 		newImage.src = base64EncString;
		// 		console.log("post newImage.src");
		// 	}) ;

		
		// console.log("pre awaitfunction()");
		// await awaitfunction;
		// console.log("post awaitfunction()");

		// console.log("pre  return m_spmimage");
		// return m_spmimage;

		//Try

		console.log("static async getSPMImageFromBitmapFile(base64EncString)");
		var m_spmimage=null;

		let newImage = new Image();
		newImage.src = base64EncString;

		let pr0 = new Promise( resolve => {
			newImage.onload = function(){

				if (newImage.width >0 && newImage.height>0){
					console.log("newImage.width >0");
					let m_spmimage0 = new SPMImage();
					let newCanvas= document.createElement('canvas');
			
					newCanvas.width = newImage.width;
					newCanvas.height = newImage.height;
					let ctx = newCanvas.getContext('2d');
					ctx.drawImage(newImage, 0, 0);
			
					//No need to fill .imageBMP with data
					m_spmimage0.xpixels = newImage.width;
					m_spmimage0.ypixels = newImage.height;
					m_spmimage0.xsize_nm=10 ;
					m_spmimage0.ysize_nm=10;
					m_spmimage0.isBitmap=true;
			
					m_spmimage0.imageAsCanvas = newCanvas;

					//return m_spmimage;
					resolve(m_spmimage0);
				}
			}
		} );

		m_spmimage= await pr0;
		
		console.log("pre  return m_spmimage");
		return m_spmimage;

	}

	static async promiseSPMImageFromBitmapFile(base64EncString){
		return new Promise( (resolve) =>{
			console.log("getSPMImageFromBitmapFilePromise(base64EncString)");
			var m_spmimage=null;

			let newImage = new Image();

			newImage.onload = function(){ //the onload function will run when setting the .src value

				if (newImage.width >0 && newImage.height>0){
					console.log("newImage.width >0");
					let m_spmimage0 = new SPMImage();
					let newCanvas= document.createElement('canvas');
			
					newCanvas.width = newImage.width;
					newCanvas.height = newImage.height;
					let ctx = newCanvas.getContext('2d');
					ctx.drawImage(newImage, 0, 0);
			
					//No need to fill .imageBMP with data
					m_spmimage0.xpixels = newImage.width;
					m_spmimage0.ypixels = newImage.height;
					m_spmimage0.xsize_nm=10 ;
					m_spmimage0.ysize_nm=10;
					m_spmimage0.isBitmap=true;
			
					m_spmimage0.imageAsCanvas = newCanvas;

					//return m_spmimage;
					resolve(m_spmimage0);
				}
			};
			
			newImage.src = base64EncString;

			
		});
	}
	
}