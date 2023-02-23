using System.Diagnostics;
using MathNet.Numerics;

var shell = SpaceDivisor<int>.HypercubeShell(6,15);

var watch = new Stopwatch();
watch.Restart();
foreach(var v in shell){
    var sum = v.Span[0]; //idk
}
System.Console.WriteLine("Time to create "+watch.ElapsedMilliseconds);