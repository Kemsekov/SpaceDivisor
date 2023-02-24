using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra.Single;
namespace SpaceDivisor.Tests;

record TestData(int id, Vector position);

public class SpaceDivisorTests
{

    public SpaceDivisorTests(){
        Dimensions = 3;
        TestData = 
            Enumerable.Range(0,1000)
            .Select(x=>
                new TestData(
                    x,
                    DenseVector.Create(Dimensions,y=>Random.Shared.NextSingle()*2-1.5f))
                ).ToList();
        SpaceDivisor = new SpaceDivisor<TestData>(TestData,x=>x.position);
        
    }

    public int Dimensions { get; }
    internal IEnumerable<TestData> TestData { get; }
    internal SpaceDivisor<TestData> SpaceDivisor { get; }

    [Fact]
    public void RightRowSize()
    {
        Assert.Equal(SpaceDivisor.RowSize,Math.Ceiling(Math.Pow(TestData.Count(),1.0f/Dimensions)));
    }
    [Fact]
    public void RightDimensions(){
        Assert.Equal(SpaceDivisor.Dimensions,Dimensions);
    }
    [Fact]
    public void RightSize()
    {
        Assert.Equal(SpaceDivisor.Size,TestData.Count());
    }
    [Fact]
    public void GetNormalized_Works(){
        var maxVector = TestData.Aggregate((n1,n2)=>{
            return new TestData(-1,(Vector)n1.position.PointwiseMaximum(n2.position));
        }).position;
        var minVector = TestData.Aggregate((n1,n2)=>{
            return new TestData(-1,(Vector)n1.position.PointwiseMinimum(n2.position));
        }).position;
        var scalar = maxVector-minVector;
        foreach(var d in TestData){
            var normalized = SpaceDivisor.GetNormalized(d.position);
            var expected = (d.position-minVector).PointwiseDivide(scalar);
            Assert.Equal(expected,normalized);
        }
    }
    [Fact]
    public void GetIndex_Works(){
        foreach(var d in TestData){
            var normalized = SpaceDivisor.GetNormalized(d.position);
            var index = new int[Dimensions];
            for(int i = 0;i<Dimensions;i++)
                index[i] = (int)MathF.Min(normalized[i]*SpaceDivisor.RowSize,SpaceDivisor.RowSize-1);
            
            var real = SpaceDivisor.GetIndex(d.position);
            Assert.Equal(index,real);
        }
    }
    [Fact]
    public void ReadByIndex_Works(){
        foreach(var d in TestData){
            var index = SpaceDivisor.GetIndex(d.position);
            var cellByIndex = SpaceDivisor[index];
            var cellByValue = SpaceDivisor[d];
            Assert.Equal(cellByIndex,cellByValue);
            Assert.Contains(d,cellByIndex);
        }
    }
    [Fact]
    public void CellContainsRightElements(){
        var grouped = TestData.GroupBy(x=>SpaceDivisor.GetIndex(x.position),new ArrayComparer<int>());
        foreach(var d in grouped){
            var first = d.First();
            var cell = SpaceDivisor[first];
            Assert.Equal(d,cell);
        }
    }
    
    [Fact]
    public void Near_Works(){
        foreach(var d in TestData)
        {
            var near = SpaceDivisor.Near(d);
            CheckNear(near);
        }
    }

    void CheckNear(IEnumerable<TestData> near)
    {
        var nearCoordinates = near
                        .Select(x => SpaceDivisor.GetIndex(x.position))
                        .Distinct(new ArrayComparer<int>());

        var expected = nearCoordinates.SelectMany(x => SpaceDivisor[x]);
        var expectedIds = expected.Select(x => x.id).Order();
        var nearIds = near.Select(x => x.id).Order();

        Assert.Equal(expectedIds, nearIds);
    }

    [Fact]
    public void ExpandedNear_Works(){
        foreach(var d in TestData){
            var near = SpaceDivisor.ExpandedNear(d,Random.Shared.Next(5));
            CheckNear(near);
        }
    }
}