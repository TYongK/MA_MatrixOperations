    static dynamic MatrixSubtract(double[][] matrixA, double[][] matrixB)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;
        int bRows = matrixB.Length; int bCols = matrixB[0].Length;
        if (aRows != bRows || aCols != bCols)
        {
            return "Non-conformable matrices";
        }

        double[][] result = MatrixCreate(aRows, bCols);

        for (int i = 0; i < aRows; ++i)
        for (int j = 0; j < bCols; ++j)
            result[i][j] = matrixA[i][j] - matrixB[i][j];

        return result;
    }

    static dynamic MatrixAddition(double[][] matrixA, double[][] matrixB)
    {
        int aRows = matrixA.Length; int aCols = matrixA[0].Length;
        int bRows = matrixB.Length; int bCols = matrixB[0].Length;
        if (aRows != bRows || aCols != bCols)
        {
            return "Non-conformable matrices";
        }

        double[][] result = MatrixCreate(aRows, bCols);

        for (int i = 0; i < aRows; ++i)
        for (int j = 0; j < bCols; ++j)
            result[i][j] = matrixA[i][j] + matrixB[i][j];

        return result;
    }
