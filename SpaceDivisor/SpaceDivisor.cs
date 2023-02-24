
using System;
using System.Buffers;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Single;
namespace MathNet.Numerics
{
    /// <summary>
    /// Divides data in some n-dimensional space into chunks and allows to look into near chunks.<br/>
    /// Can be used to find close elements in data in constant time.<br/>
    /// Works well for cases where number of dimensions is small, but number of elements in data is big
    /// and data is uniformly distributed.<br/>
    /// </summary>
    /// <typeparam name="TValue"></typeparam>
    public partial class SpaceDivisor<TValue> : IDisposable
    {
        /// <summary>
        /// Creates new instance of space divisor from given data and function to get vector from particular element<br/>
        /// Creation takes <see langword="O(n)"/> time.
        /// </summary>
        /// <param name="data"></param>
        /// <param name="getPosition"></param>
        public SpaceDivisor(IEnumerable<TValue> data, Func<TValue, Vector> getPosition)
        {
            var size = data.Count();
            Data = data;
            GetPosition = getPosition;
            Dimensions = getPosition(data.First()).Count;
            RowSize = (int)MathF.Ceiling(MathF.Pow(size, 1.0f / Dimensions));
            var toAllocate = (int)Math.Pow(RowSize, Dimensions);
            //we allocate in roughly O(n) space complexity
            storage = ArrayPool<List<TValue>>.Shared.Rent(toAllocate);
            normalizer = Normalizer();
            FillStorage();
        }
        /// <summary>
        /// Computes index of a cell where given vector should be
        /// </summary>
        public int[] GetIndex(Vector v)
        {
            var index = GetNormalized(v);
            var multiplier = Math.Max(RowSize, 1);
            var result = new int[Dimensions];
            for (int i = 0; i < Dimensions; i++)
            {
                result[i] = (int)MathF.Min(index[i] * multiplier,RowSize-1);
            }
            return result;
        }
        /// <summary>
        /// Does not create new list when reading not initialized cell
        /// </summary>
        List<TValue>? FastRead(int[] index)
        {
            var intIndex = 0;
            for (int i = 0; i < Dimensions; i++)
            {
                if (index[i] >= RowSize || index[i] < 0) return null;
                intIndex += (int)(index[i] * MathF.Pow(RowSize, i));
            }
            return storage[intIndex];
        }
        /// <summary>
        /// Returns cell under given index
        /// </summary>
        public IEnumerable<TValue> this[params int[] index]
        {
            get
            {
                if (index.Length != Dimensions)
                    throw new ArgumentException($"Index have not enough dimensions! Given {index.Length}. Required {Dimensions}");
                var intIndex = 0;
                for (int i = 0; i < Dimensions; i++)
                {
                    if (index[i] >= RowSize || index[i] < 0)
                    {
                        throw new IndexOutOfRangeException($"{index[i]}>={RowSize} or {index[i]}<0");
                    }
                    intIndex += (int)(index[i] * MathF.Pow(RowSize, i));
                }
                var value = storage[intIndex];
                if (value is null)
                {
                    return Enumerable.Empty<TValue>();
                }
                return value;
            }
        }
        List<TValue> GetOrCreate(int[] index)
        {
                if (index.Length != Dimensions)
                    throw new ArgumentException($"Index have not enough dimensions! Given {index.Length}. Required {Dimensions}");
                var intIndex = 0;
                for (int i = 0; i < Dimensions; i++)
                {
                    if (index[i] >= RowSize || index[i] < 0)
                    {
                        throw new IndexOutOfRangeException($"{index[i]}>={RowSize} or {index[i]}<0");
                    }
                    intIndex += (int)(index[i] * MathF.Pow(RowSize, i));
                }
                var value = storage[intIndex];
                if (value is null)
                {
                    value = new List<TValue>();
                    storage[intIndex] = value;
                }
                return value;
        }
        /// <summary>
        /// Combines elements from cell where given element is present and elements from near cells<br/>
        /// Two cells are near if they differ only by one coordinate by <see langword="1"/>.<br/>
        /// For two dimensional case it will merge following cells:<br/>
        /// If given value is under cell <see langword="[2,2]"/>, then <br/>
        /// return all values from cells <see langword="[2,2],[1,2],[3,2],[2,1],[2,3]"/> (like cross)<br/>
        /// For three dimensional case:<br/>
        /// If given value is under cell <see langword="[2,2,2]"/>, then <br/>
        /// return all values from cells <see langword="[2,2,2],[1,2,2],[3,2,2],[2,1,2],[2,3,2],[2,2,1],[2,2,3]"/>(like 3-dim cross)
        /// </summary>
        public IEnumerable<TValue> Near(TValue v) => Near(GetPosition(v));
        ///<inheritdoc cref="Near(TValue)"/>
        public IEnumerable<TValue> Near(Vector v)
        {
            var index = GetIndex(v);

            foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;

            for (int i = 0; i < Dimensions; i++)
            {
                var originalValue = index[i];
                index[i] = originalValue + 1;

                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;

                index[i] = originalValue - 1;
                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;
                
                index[i] = originalValue;
            }
        }
        /// <summary>
        /// Enumerates N-dimensional space.<br/>
        /// It means that for number in base <see langword="k"/> of length <see langword="dim"/><br/>
        /// It will return all values from zero to maximum.<br/>
        /// For example with <see langword="dim=2, k=10"/> it will return 
        /// all numbers from <see langword="00"/> up to <see langword="99"/>
        /// </summary>
        public static IEnumerable<ReadOnlyMemory<int>> EnumerateNDimensionalSpace(int dim, int k)
        {
            var index = new int[dim];
            yield return index;
            while (true)
            {
                index[0] += 1;
                for (int i = 0; i < dim; i++)
                {
                    if (index[i] >= k)
                    {
                        if (i + 1 >= dim) yield break;
                        index[i + 1] += 1;
                        index[i] -= k;
                    }
                    else break;
                }
                yield return index;
            }
        }
        /// <summary>
        /// For <see langword="dim"/>-dimensional hypercube, divide it into <see langword="k^dim"/>
        /// chunks and return only those chunks that touches surface of hypercube.<br/>
        /// Example from <see langword="dim=2, k=4"/><br/>
        /// #--#--#--# ---> #--#--#--# <br/>
        /// #--#--#--# ---> #----------# <br/>
        /// #--#--#--# ---> #----------# <br/>
        /// #--#--#--# ---> #--#--#--# <br/>
        /// I hope you smart enough to expand this idea to higher dimensions
        /// </summary>
        public static IEnumerable<ReadOnlyMemory<int>> HypercubeShell(int dim, int k)
        {
            var index = new int[dim];
            if (dim == 1)
            {
                yield return index;
                index[0] = k - 1;
                yield return index;
                yield break;
            }

            var cubeShellSlice = HypercubeShell(dim - 1, k);
            var hyperplane = EnumerateNDimensionalSpace(dim - 1, k);
            foreach (var v in hyperplane)
            {
                v.Span.CopyTo(index);
                yield return index;
            }

            Array.Fill(index, 0);
            foreach (var v in cubeShellSlice)
            {
                v.Span.CopyTo(index);
                for (int i = 1; i < k - 1; i++)
                {
                    index[^1] = i;
                    yield return index;
                }
            }

            index[^1] = k - 1;
            foreach (var v in hyperplane)
            {
                v.Span.CopyTo(index);
                yield return index;
            }
        }
        
        ///<inheritdoc cref="ExpandedNear(Vector,int)"/>
        public IEnumerable<TValue> ExpandedNear(TValue v, int radius = 1) => ExpandedNear(GetPosition(v), radius);
        /// <summary>
        /// Allows you to take layers of near cells individually.<br/>
        /// For example if you searching for closest element to given value, start from <see langword="radius=0"/>(aka cell where element is stored),
        /// and if there is no elements except value itself, increase radius and you will get elements from near cells,
        /// but not from the previous ones.<br/>
        /// So you can 'expand' your search depending of whether you found what you need or not.<br/>
        /// It is guaranteed that return of this function with higher <see langword="radius"/> will result
        /// in elements that is farther than any elements from smaller <see langword="radius"/><br/>
        /// I also recommend you to sort output of this function, because it gives results from chunks,
        /// but in same chunks order of closest elements to given value is not defined.<br/>
        /// For <see langword="radius=0"/><br/>
        /// --------<br/>
        /// ---#---<br/>
        /// --------<br/>
        /// For <see langword="radius=1"/><br/>
        /// #--#--#<br/>
        /// #------#<br/>
        /// #--#--#<br/>
        /// For <see langword="radius=2"/><br/>
        /// #--#--#--#<br/>
        /// #----------#<br/>
        /// #----------#<br/>
        /// #--#--#--#<br/>
        /// And so on, but for any n-dimensional cube, in general returns closest cells in hypercube shell shape
        /// </summary>
        public IEnumerable<TValue> ExpandedNear(Vector v, int radius = 1)
        {
            var index = GetIndex(v);
            if (radius == 0)
                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;

            if (radius <= 0) yield break;

            var center = radius + 1;
            foreach (var position in HypercubeShell(Dimensions, 2 * radius + 1))
            {
                var pos = position.Span;
                for (int i = 0; i < pos.Length; i++)
                {
                    index[i] += pos[i] - center + 1;
                }
                foreach (var value in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return value;
                for (int i = 0; i < position.Span.Length; i++)
                {
                    index[i] -= position.Span[i] - center + 1;
                }
            }
        }
        /// <returns>Values stored in same cell as given value</returns>
        public IEnumerable<TValue> this[TValue value] => this[GetPosition(value)];
        ///<inheritdoc cref="this[TValue]"/>
        public IEnumerable<TValue> this[Vector value] => this[GetIndex(value)];
        void FillStorage()
        {
            foreach (var v in Data)
            {
                GetOrCreate(GetIndex(GetPosition(v))).Add(v);
            }
        }

        (Vector scalar, Vector minVector) Normalizer()
        {
            var m1 = new float[Dimensions];
            Vector maxVector = new DenseVector(m1);
            Array.Fill(m1, float.MinValue);
            var m2 = new float[Dimensions];
            Vector minVector = new DenseVector(m2);
            Array.Fill(m2, float.MaxValue);

            foreach (var v in Data)
            {
                maxVector = (Vector)maxVector.PointwiseMaximum(GetPosition(v));
                minVector = (Vector)minVector.PointwiseMinimum(GetPosition(v));
            }
            var scalar = maxVector - minVector;
            for (int i = 0; i < Dimensions; i++)
                //to avoid division by zero
                if (Math.Abs(scalar[i]) <= float.Epsilon)
                    scalar[i] = 1;
            return ((Vector)scalar, minVector);
        }
        /// <returns>Vector normalized in respect to given to current space divisor data</returns>
        public Vector GetNormalized(Vector v)
        {
            return (Vector)(v - normalizer.minVector).PointwiseDivide(normalizer.scalar);
        }
        /// <summary>
        /// When destroyed also disposes
        /// </summary>
        ~SpaceDivisor()
        {
            Dispose();
        }
        bool disposed = false;
        ///<inheritdoc/>
        public void Dispose()
        {
            lock (storage)
            {
                if (disposed) return;
                ArrayPool<List<TValue>>.Shared.Return(storage, true);
                disposed = true;
            }
        }
        /// <summary>
        /// Number of dimensions
        /// </summary>
        public readonly int Dimensions;
        /// <summary>
        /// Row size of each dimension
        /// </summary>
        public readonly int RowSize;
        private readonly List<TValue>[] storage;
        private (Vector scalar, Vector minVector) normalizer;
        /// <summary>
        /// Data that populate current space divisor
        /// </summary>
        public IEnumerable<TValue> Data { get; }
        /// <summary>
        /// Function to get position of data elements with respect to current space dimensions
        /// </summary>
        public Func<TValue, Vector> GetPosition { get; }
    }
}