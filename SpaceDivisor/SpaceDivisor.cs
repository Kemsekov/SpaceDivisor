
using System;
using System.Buffers;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Single;
namespace MathNet.Numerics
{
    public partial class SpaceDivisor<TValue> : IDisposable
    {
        public SpaceDivisor(IEnumerable<TValue> data, Func<TValue, Vector> getPosition)
        {
            Size = data.Count();
            Data = data;
            GetPosition = getPosition;
            Dimensions = getPosition(data.First()).Count;
            RowSize = (int)MathF.Ceiling(MathF.Pow(Size, 1.0f / Dimensions));
            var toAllocate = (int)Math.Pow(RowSize, Dimensions);
            //we allocate in roughly O(n) space complexity
            storage = ArrayPool<List<TValue>>.Shared.Rent(toAllocate);
            normalizer = Normalizer();
            FillStorage();
        }
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
        public IEnumerable<TValue> Near(TValue v) => Near(GetPosition(v));
        public IEnumerable<TValue> Near(Vector v)
        {
            var index = GetIndex(v);

            foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                if (GetPosition(v1) != v)
                    yield return v1;

            for (int i = 0; i < Dimensions; i++)
            {
                var originalValue = index[i];

                index[i] = originalValue + 1;
                index[i] = Math.Min(index[i], RowSize - 1);
                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;
                index[i] = originalValue - 1;
                index[i] = Math.Max(index[i], 0);
                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    yield return v1;
                index[i] = originalValue;
            }
        }
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
        public IEnumerable<TValue> ExpandedNear(TValue v, int radius = 1) => ExpandedNear(GetPosition(v), radius);
        public IEnumerable<TValue> ExpandedNear(Vector v, int radius = 1)
        {
            var index = GetIndex(v);
            if (radius == 0)
                foreach (var v1 in FastRead(index) ?? Enumerable.Empty<TValue>())
                    if (GetPosition(v1) != v)
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

        public IEnumerable<TValue> this[TValue index] => this[GetPosition(index)];
        public IEnumerable<TValue> this[Vector index] => this[GetIndex(index)];
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
        public Vector GetNormalized(Vector v)
        {
            return (Vector)(v - normalizer.minVector).PointwiseDivide(normalizer.scalar);
        }
        ~SpaceDivisor()
        {
            Dispose();
        }
        bool disposed = false;
        public void Dispose()
        {
            lock (storage)
            {
                if (disposed) return;
                ArrayPool<List<TValue>>.Shared.Return(storage, true);
                disposed = true;
            }
        }

        public readonly int Dimensions;
        public readonly int RowSize;
        public readonly int Size;
        private readonly List<TValue>[] storage;
        private (Vector scalar, Vector minVector) normalizer;

        public IEnumerable<TValue> Data { get; }
        public Func<TValue, Vector> GetPosition { get; }
    }
}