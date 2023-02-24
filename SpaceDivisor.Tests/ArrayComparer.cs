class ArrayComparer<T> : IEqualityComparer<T[]>
{
    public bool Equals(T[]? x, T[]? y)
    {
        if (x is null || y is null || x.Length != y.Length)
        {
            return false;
        }

        for (int i = 0; i < x.Length; i++)
        {
            if (!EqualityComparer<T>.Default.Equals(x[i], y[i]))
            {
                return false;
            }
        }

        return true;
    }

    public int GetHashCode(T[] obj)
    {
        var h = new HashCode();
        foreach(var v in obj)
            h.Add(v);

        return h.ToHashCode();
    }
}