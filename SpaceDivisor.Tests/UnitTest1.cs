using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Single;
namespace SpaceDivisor.Tests;

record TestData(int id, Vector position);

public class UnitTest1
{

    [Fact]
    public void Test1()
    {
        var testData = 
            Enumerable.Range(0,1000)
            .Select(x=>
                new TestData(
                    x,
                    DenseVector.Create(2,y=>Random.Shared.NextSingle()))
                );

        var d = new SpaceDivisor<TestData>(testData,x=>x.position);
        
    }
}