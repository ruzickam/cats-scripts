#!/path/to/cats

// load topology
var mp = new Topology("dna_joined.parm7");
mp.printInfo();

// init trajectory
var mt = new TrajPool(mp)
mt.addTrajListFrom(1,20,"../client01/storage");
mt.addTrajListFrom(1,20,"../client02/storage");
mt.addTrajListFrom(1,20,"../client03/storage");
mt.addTrajListFrom(1,20,"../client04/storage");
mt.addTrajListFrom(1,20,"../client05/storage");
mt.printInfo();

// subsystem
var sl = new Selection(mp,":1-30");
// SET BP1
var bp1residA = 7;
var bp1residB = 22;

// helper objects
var ms = new Snapshot(mp);

var of = new OFile("dat.log");
of.printf("%5d %14.3f %14.3f %14.3f %14.3f %14.3f\n","snap","sh/pmflib","sh/3dna","rmsd","d1/O6-N3","d2/N1-O2");

// init 3dna and PMFLib
PMFLib.init(ms,"cvs.in");
var x3dna = new x3DNA();
x3dna.setParameterType("local");
var nastat = new NAStat();

printf("\n# analyzing ...\n");

// analyze trajectory
var i = 0;
while( mt.read(ms) ){
   i++;
   mt.printProgress();
   x3dna.analyze(ms,sl);
   PMFLib.setCoordinates(ms);
   nastat.addSample(x3dna);
   var a1 = PMFLib.getCVValue(0);
   var s1 = PMFLib.getCVValue(1);
   var d1 = PMFLib.getCVValue(2);
   var d2 = PMFLib.getCVValue(3);
   var bp1index = x3dna.getBPIndex(bp1residA,bp1residB);
   if( bp1index == -1 ){
    of.printf("%5d %14.3f %14.3f %14.3f %14.3f %14.3f\n",i,s1,"?",a1,d1,d2);
   } else {
    if( x3dna.areBPParamsValid(bp1index) == true ){
     var s2 = x3dna.getBPShear(bp1index);
     of.printf("%5d %14.3f %14.3f %14.3f %14.3f %14.3f\n",i,s1,s2,a1,d1,d2);
    } else {
     of.printf("%5d %14.3f %14.3f %14.3f %14.3f %14.3f\n",i,s1,"-",a1,d1,d2);
    }
   }
}

of.close();
nastat.printResults("na.stat");
