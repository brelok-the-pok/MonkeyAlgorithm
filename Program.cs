using System;
using System.IO;
using static System.Math;

namespace MonkeyAlgorithm
{
    class Program
    {
        private static Random random = new Random();
        private const int rank = 2;

        static void Main()
        {

            double[] min = new double[rank] { -10, 1};//Минимум функции для проверки 
            Function f = new Function(BookinN6, 2, min);

            double[][] populationPosition;//Текущие позиции каждой особи
            double[] globalMin = new double[rank];//Глобальный минимум
            double[] gravity;//Центр тяжести функции
            double[][] interval = new double[rank][];//Интервал существования функции
            
            interval[0] = new double[2] { -15, -5 };
            interval[1] = new double[2] { -3, 3 };

            double eps = 0.0001;//Некоторая точность
            double b = 2;//1;//Длина локального прыжка
            double v = 4;//2;//Длина глобального прыжка

            int S = 50;//Размер популяции
            int iteration = 0;//Номер итерации
            int iterationsMax = 1000;//Максимальное число итераций

            populationPosition = InitializePopulation(interval, S);//Инициализация популяции обезьян

            globalMin = (double[])populationPosition[0].Clone();

            while (iteration < iterationsMax)
            {
                iteration++;

                for (int i = 0; i < S; i++)
                {
                    gravity = CalculateGravity(populationPosition, S);//Нахождение центра тяжести функции 

                    populationPosition[i] = MethodNR(f, populationPosition[i], eps, interval[0]);//Процесс движения вверх

                    globalMin = CheckGlobalMin(f, populationPosition, globalMin, S);//Сохранение глобального минимума

                    //Локальные прыжки
                    populationPosition[i] = WatchJumpProcess(f, populationPosition[i], i, interval, b);

                    //Глобальные прыжки
                    populationPosition[i] = SomeResaultProcess(populationPosition[i], interval, gravity, v);

                    if (GlobalIsMinFound(f, globalMin, f.minimum, eps))
                    {
                        WriteABeautifullAnswer($"На итерация номер {iteration} был найдем минимум функции особью {i}", f[globalMin], globalMin);
                        return;
                    }
                }
            }

            Logger.WriteLine($"Число итераций исчерпано. {S} особей ни смогли найти минимум за {iteration} итераций\n, точный минимум функции с точность {eps} не найден, последние глобальные значения:\n" +
                           $"Значение функции {f[globalMin]} при x = {globalMin[0]} и y = {globalMin[1]}");
            Logger.WriteLine($"Реальное значение минимума функции:\n" +
                           $"Значение функции {f[min]} при x = {min[0]} и y = {min[1]}");
            WriteABeautifullAnswer($"Глобальный минимум не был найден", f[globalMin], globalMin);

        }
        private static bool GlobalIsMinFound(Function f, double[] min, double[] realMin, double eps)
        {
            if(Abs(f[min] - f[realMin]) <= eps)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        private static double[][] InitializePopulation(double[][] interval, int S)
        {
            double[][] population = new double[S][];
            for (int i = 0; i < S; i++)
            {
                population[i] = new double[rank];
            }
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < S; j++)
                {
                    population[j][i] = random.NextDouble() * (interval[i][1] - interval[i][0]) + interval[i][0];
                }
            }
            return population;
        }
        private static double[] CalculateGravity(double[][] currentPosition, int S)
        {
            double[] gravity = new double[rank];
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < S; j++)
                {
                    gravity[i] += currentPosition[j][i];
                }
                gravity[i] /= S;
            }
            return gravity;
        }
        private static double[] WatchJumpProcess(Function f, double[] currentPosition,
            int individual, double[][] interval, double jumpLength)
        {
            double[] nextPosition = (double[])currentPosition.Clone();

            double randomValue;
            for (int i = 0; i < rank; i++)
            {
                do
                {
                    randomValue = random.NextDouble() * 2 * jumpLength + (currentPosition[i] - jumpLength);
                } while (randomValue < interval[i][0] || randomValue > interval[i][1]);

                nextPosition[i] = randomValue;
                if (f[nextPosition] > f[currentPosition])
                {
                    i--;
                }
            }
            return nextPosition;
        }
        private static double[] SomeResaultProcess(double[] currentPosition, double[][] interval, double[] gravity, double jumpLength)
        {
            double[] nextPosition = (double[])currentPosition.Clone();
            double randomValue;
            double x;

            for (int i = 0; i < rank; i++)
            {
                do
                {
                    randomValue = random.NextDouble() * 2 * jumpLength - jumpLength;
                    x = currentPosition[i] + randomValue * (gravity[i] - currentPosition[i]);
                } while (x < interval[i][0] || x > interval[i][1]);

                nextPosition[i] = x;
            }
            return nextPosition;
        }
        private static double[] CheckGlobalMin(Function f, double[][] currentPosition, double[] currentGlobalMin, int S)
        {
            for (int i = 0; i < S; i++)
            {
                if (f[currentGlobalMin] > f[currentPosition[i]])
                {
                    currentGlobalMin = (double[])currentPosition[i].Clone();
                }
            }
            return currentGlobalMin;
        }
        static double BookinN6(params double[] x)//x[-15;-5]  y[-3;3]x = {-10, 1} f = 0
        {

            return 100 * Sqrt(Abs(x[1] - 0.01 * x[0] * x[0]) ) + 0.01* Abs(x[0] + 10);
        }
        static double CrossInTrayFunction(params double[] x)
        {
            double exp = Exp(Abs( 100 - (Sqrt(x[0] * x[0] + x[1] * x[1]) / PI)));
            double sinsin = Sin(x[0]) * Sin(x[1]);

            return -0.0001 * Pow(Abs(sinsin * exp) + 1, 0.1);
        }
        static double ShafferFunction4(params double[] x)
        {
            double num1 = Pow(Cos(Sin(Abs(x[0] * x[0] - x[1] * x[1]))), 2) - 0.5;
            double num2 = Pow(1 + 0.001 * (x[0] * x[0] + x[1] * x[1]), 2);

            return 0.5 + num1 / num2;
        }

        static double EngholdFunc(params double[] x)
        {
            return - ( x[1] + 47) * Sin(Sqrt(Abs(x[0] / 2 + (x[1] + 47)))) - x[0] * Sin(Sqrt(Abs(x[0] - (x[1] + 47))));
        }

        static public double[] MethodNR(Function f, double[] startPos, double eps, double[] interval)//Метод Ньютона-Рафсона
        {
            int k = 0;
            int rank = f.Rank;
            double[] gradVector;
            double[] curX = (double[])startPos.Clone();
            double[,] hessian;//Матрица Гессе
            double lymbda;
            while (k < 100)
            {
                gradVector = Gradient(f, curX);
                if (EuclideanNorm(gradVector) <= eps)
                {
                    return curX;
                }

                hessian = FillHessian(f, curX);
                double[,] inverceHessian = MatrixInverse(hessian);
                for (int i = 0; i < rank; i++)
                {
                    for (int j = 0; j < rank; j++)
                    {
                        inverceHessian[i, j] = -inverceHessian[i, j];
                    }
                }

                double[,] searchDirection;//Направление поиска

                double[,] gradForMatixMul = new double[rank, rank];
                for (int i = 0; i < rank; i++)
                {
                    gradForMatixMul[i, 0] = gradVector[i];
                }

                searchDirection = MatrixMul(inverceHessian, gradForMatixMul);

                for (int i = 0; i < rank; i++)
                {
                    if (searchDirection[i, 0] != 0)
                    {
                        double[] z;
                        lymbda = 2;
                        do
                        {
                            lymbda /= 2;
                            z = (double[])curX.Clone();
                            z[i] += lymbda * searchDirection[i, 0];
                        } while (Abs(f[z] - f[curX]) > eps && (f[z] > f[curX]));

                        if(notInInterval(z, interval))
                        {
                            return curX;
                        }

                        curX[i] += lymbda * searchDirection[i, 0];
                    }
                }
                k++;
            }
            return curX;
        }
        public static bool notInInterval(double[] x, double[] interval)
        {
            foreach (double xval in x)
            {
                if(xval > interval[1] || xval < interval[0])
                {
                    return true;
                }
            }
            return false;
        }
        public static double[,] FillHessian(Function f, double[] x)
        {
            int rank = f.Rank;
            double[,] hessian = new double[rank, rank];
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < rank; j++)
                {
                    hessian[i, j] = secondPatrialDerivative(f, x, i, j);
                }
            }
            return hessian;
        }
        public static double[,] MatrixMul(double[,] matrix1, double[,] matrix2)//Произведение двух матриц
        {
            if (matrix1.GetLength(0) != matrix2.GetLength(1))
            {
                Console.WriteLine("Неверные размерности матриц");
                return new double[matrix1.GetLength(0), matrix1.GetLength(1)];
            }
            double[,] resMatrix = new double[matrix1.GetLength(1), matrix2.GetLength(0)];
            for (int i = 0; i < matrix1.GetLength(0); i++)
            {
                for (int j = 0; j < matrix2.GetLength(1); j++)
                {
                    for (int k = 0; k < matrix2.GetLength(0); k++)
                    {
                        resMatrix[i, j] += matrix1[i, k] * matrix2[k, j];
                    }
                }
            }
            return resMatrix;
        }
        public static double[,] MatrixInverse(double[,] matrix)//Нахождение обратно мартицы
        {
            int rank = matrix.GetLength(0);
            double[,] reversedMatrix = (double[,])matrix.Clone();
            double det = MatrixDeterminant(reversedMatrix);
            double[,] minorMarix = new double[rank, rank];
            for (int i = 0; i < rank; i++)
            {
                for (int j = 0; j < rank; j++)
                {
                    double[,] minor = new double[rank - 1, rank - 1];
                    for (int a = 0, a1 = 0; a < rank; a++)
                    {
                        int b1 = 0;
                        for (int b = 0; b < rank; b++)
                        {
                            if (a != i && b != j)
                            {
                                minor[a1, b1] = reversedMatrix[a, b];
                                b1++;
                            }
                        }
                        if (a != i)
                            a1++;
                    }
                    minorMarix[i, j] = Pow(-1.0, i + j) * MatrixDeterminant(minor);
                }
            }
            reversedMatrix = MatrixTranspose(minorMarix);
            reversedMatrix = MatrixOnConst(reversedMatrix, 1.0 / det);

            return reversedMatrix;
        }
        public static double[] Gradient(Function f, double[] x)
        {
            double delta = 1e-8;
            int rank = f.Rank;
            double[] gradArr = new double[rank];
            double[] curX = (double[])x.Clone();
            for (int i = 0; i < rank; i++)
            {
                curX[i] = x[i] + delta;
                gradArr[i] = (f[curX] - f[x]) / delta;
                curX[i] = x[i];
            }
            return gradArr;
        }
        public static double EuclideanNorm(params double[] vals)
        {
            double norm = 0;
            foreach (double val in vals) norm += Pow(val, 2);
            return Sqrt(norm);
        }
        public static double MatrixDeterminant(double[,] matrix)//Нахождение определителя матрицы 
        {

            if (matrix == null)
            {
                throw new ArgumentNullException("MatrixDeterminant(double[,] matrix) = null");
            }
            double det;
            if (matrix.GetLength(0) == 1)
            {
                det = matrix[0, 0];
            }
            else if (matrix.GetLength(0) == 2)
            {
                det = (matrix[0, 0] * matrix[1, 1]) - (matrix[0, 1] * matrix[1, 0]);
            }
            else if (matrix.GetLength(0) == 3)
            {
                det =
                matrix[0, 0] * matrix[1, 1] * matrix[2, 2] +
                matrix[0, 1] * matrix[1, 2] * matrix[2, 0] +
                matrix[0, 2] * matrix[1, 0] * matrix[2, 1] -
                matrix[0, 2] * matrix[1, 1] * matrix[2, 0] -
                matrix[0, 0] * matrix[1, 2] * matrix[2, 1] -
                matrix[0, 1] * matrix[1, 0] * matrix[2, 2];

            }
            else if (matrix.GetLength(0) == 4)
            {
                det =
                matrix[0, 3] * matrix[1, 2] * matrix[2, 1] * matrix[3, 0] - matrix[0, 2] * matrix[1, 3] * matrix[2, 1] * matrix[3, 0] -
                matrix[0, 3] * matrix[1, 1] * matrix[2, 2] * matrix[3, 0] + matrix[0, 1] * matrix[1, 3] * matrix[2, 2] * matrix[3, 0] +
                matrix[0, 2] * matrix[1, 1] * matrix[2, 3] * matrix[3, 0] - matrix[0, 1] * matrix[1, 2] * matrix[2, 3] * matrix[3, 0] -
                matrix[0, 3] * matrix[1, 2] * matrix[2, 0] * matrix[3, 1] + matrix[0, 2] * matrix[1, 3] * matrix[2, 0] * matrix[3, 1] +
                matrix[0, 3] * matrix[1, 0] * matrix[2, 2] * matrix[3, 1] - matrix[0, 0] * matrix[1, 3] * matrix[2, 2] * matrix[3, 1] -
                matrix[0, 2] * matrix[1, 0] * matrix[2, 3] * matrix[3, 1] + matrix[0, 0] * matrix[1, 2] * matrix[2, 3] * matrix[3, 1] +
                matrix[0, 3] * matrix[1, 1] * matrix[2, 0] * matrix[3, 2] - matrix[0, 1] * matrix[1, 3] * matrix[2, 0] * matrix[3, 2] -
                matrix[0, 3] * matrix[1, 0] * matrix[2, 1] * matrix[3, 2] + matrix[0, 0] * matrix[1, 3] * matrix[2, 1] * matrix[3, 2] +
                matrix[0, 1] * matrix[1, 0] * matrix[2, 3] * matrix[3, 2] - matrix[0, 0] * matrix[1, 1] * matrix[2, 3] * matrix[3, 2] -
                matrix[0, 2] * matrix[1, 1] * matrix[2, 0] * matrix[3, 3] + matrix[0, 1] * matrix[1, 2] * matrix[2, 0] * matrix[3, 3] +
                matrix[0, 2] * matrix[1, 0] * matrix[2, 1] * matrix[3, 3] - matrix[0, 0] * matrix[1, 2] * matrix[2, 1] * matrix[3, 3] -
                matrix[0, 1] * matrix[1, 0] * matrix[2, 2] * matrix[3, 3] + matrix[0, 0] * matrix[1, 1] * matrix[2, 2] * matrix[3, 3];
            }
            else
            {
                det = 0;
            }
            return det;
        }
        public static double[,] MatrixOnConst(double[,] matrix, double constValue)//Нахождение произведения матрицы на константу
        {
            double[,] buffer = new double[matrix.GetLength(0), matrix.GetLength(1)];
            for (int i = 0; i < buffer.GetLength(0); i++)
            {
                for (int j = 0; j < buffer.GetLength(1); j++)
                {
                    buffer[i, j] = matrix[i, j] * constValue;
                }
            }
            return buffer;
        }
        public static double[,] MatrixTranspose(double[,] matrix)//Транспонирование матрицы
        {
            double[,] buffer = (double[,])matrix.Clone();
            for (int i = 0; i < buffer.GetLength(0); i++)
            {
                for (int j = i; j < buffer.GetLength(1); j++)
                {
                    double buff = buffer[i, j];
                    buffer[i, j] = matrix[j, i];
                    buffer[j, i] = buff;
                }
            }
            return buffer;
        }
        public static double secondPatrialDerivative(Function f, double[] x, int dx1, int dx2)
        {
            double result = 0;
            double delta = 1e-4;
            double[] curX = (double[])x.Clone();
            if (dx1 == dx2)
                curX[dx1] = x[dx1] + 2 * delta;
            else
            {
                curX[dx1] = x[dx1] + delta;
                curX[dx2] = x[dx2] + delta;
            }
            double result1 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1];
            else
            {
                curX[dx1] = x[dx1] + delta;
                curX[dx2] = x[dx2] - delta;
            }
            double result2 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1];
            else
            {
                curX[dx1] = x[dx1] - delta;
                curX[dx2] = x[dx2] + delta;
            }
            double result3 = f[curX];
            if (dx1 == dx2)
                curX[dx1] = x[dx1] - 2 * delta;
            else
            {
                curX[dx2] = x[dx2] - delta;
            }
            double result4 = f[curX];
            result = (result1 - result2 - result3 + result4) / (4 * Pow(delta, 2));
            return result;
        }
        public static void WriteABeautifullAnswer(string msg, double functionValue, params double[] values)
        {
            Console.WriteLine(msg);
            string text = "";
            string text2 = "|";
            int count = 0;
            foreach (double val in values)
            {
                count += val.ToString("F4").Length;
            }
            count += functionValue.ToString("F4").Length;
            int spacesPerX = count / (values.Length + 1);
            for (int i = 0; i < count + 3 * values.Length + 6; i++)
            {
                text += "-";
            }

            Console.WriteLine(text);
            for (int i = 0; i < ((count + 3 * values.Length + 6) / 2) - (1 + "Значение функции в заданных точках".Length / 2); i++)
            {
                text2 += " ";
            }
            text2 += "Значение функции в заданных точках";
            for (int i = 0; i < ((count + 3 * values.Length + 6) / 2) - (1 + "Значение функции в заданных точках".Length / 2); i++)
            {
                text2 += " ";
            }
            text2 += "|";
            Console.WriteLine(text2);
            Console.WriteLine(text);
            string text3 = "| ";
            for (int i = 0; i < values.Length; i++)
            {
                for (int j = 0; j < (spacesPerX / 2) - 1; j++)
                {
                    text3 += " ";
                }
                text3 += string.Format("x{0}", i);
                for (int j = 0; j < (spacesPerX / 2); j++)
                {
                    text3 += " ";
                }
                text3 += "|";
            }
            for (int j = 0; j < (spacesPerX / 2) - 1; j++)
            {
                text3 += " ";
            }
            text3 += " f";
            for (int j = 0; j < (spacesPerX / 2); j++)
            {
                text3 += " ";
            }
            text3 += " |";

            Console.WriteLine(text3);
            Console.WriteLine(text);
            string text4 = "|";
            for (int i = 0; i < values.Length; i++)
            {
                for (int j = 0; j < (spacesPerX / 2) - values[i].ToString("F4").Length / 2; j++)
                {
                    text4 += " ";
                }
                text4 += values[i].ToString("F4");
                for (int j = 0; j < (spacesPerX / 2) - values[i].ToString("F4").Length / 2; j++)
                {
                    text4 += " ";
                }
                text4 += " | ";
            }
            for (int j = 0; j < (spacesPerX / 2) - functionValue.ToString("F4").Length / 2; j++)
            {
                text4 += " ";
            }
            text4 += functionValue.ToString("F4");
            for (int j = 0; j < (spacesPerX / 2) - functionValue.ToString("F4").Length / 2; j++)
            {
                text4 += " ";
            }
            text4 += " |";
            Console.WriteLine(text4);
            Console.WriteLine(text);
        }
    }

    public class Function
    {
        public delegate double function(params double[] values);
        public Function(function f, int rank, double[] minimum)
        {
            this.f = f;
            Rank = rank;
            this.minimum = minimum;
        }
        public int Rank { get; set; }//Размерность функции
        function f;
        public double[] minimum;
        virtual public double this[params double[] values]
        {
            get
            {
                return f(values);
            }
        }

    }

    public class Logger
    {
        public static void WriteLine(string text)
        {
            using (StreamWriter w = File.AppendText(Directory.GetCurrentDirectory() + @"\log.txt"))
            {
                Log(text, w);
            }
        }
        public static void Log(string logMessage, TextWriter w)
        {
            w.WriteLine($"{logMessage}");
            w.WriteLine("-------------------------------------------------------------------");
        }
        public static void DumpLog(StreamReader r)
        {
            string line;
            while ((line = r.ReadLine()) != null)
            {
                Console.WriteLine(line);
            }
        }
    }
    
}