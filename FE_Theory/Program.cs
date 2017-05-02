using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Accord.Math;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;

namespace FE_Theory
{
    class Program
    {
        static int nodeCount = 3;                      //Total # of Nodes 
        static int elementCount = 3;                   //Total # of Elements
        static int numberOfNodesPerElement = 2;        //# of Nodes/Element
        static int numberDegreesOfFreedomPerNode = 2;  //#DOF per node
        private static int activeDoF = 0;

        //XY Coordinates of the nodes
        private static Matrix<double> nodeCoordinates = Matrix<double>.Build.Dense(nodeCount, numberDegreesOfFreedomPerNode);
        //Element Connectivity Matrix
        private static Matrix<double> elementConnections = Matrix<double>.Build.Dense(elementCount, 2);
        //Element Properties
        private static Matrix<double> elementProperties = Matrix<double>.Build.Dense(elementCount, 2);

        //Matrix of Nodal Freedom (nf[nnd,nodof])
        private static Matrix<double> nodalFreedom = Matrix<double>.Build.Dense(nodeCount, numberDegreesOfFreedomPerNode);

        //Applied Loads
        private static Matrix<double> appliedLoad = Matrix<double>.Build.Dense(nodeCount, numberDegreesOfFreedomPerNode);



        static void Main(string[] args)
        {
            InitializeNodeElementMatrices();
            InitBCMatrix();
            

            for (int i = 0; i < nodalFreedom.RowCount; i++)
            {
                for (int j = 0; j < nodalFreedom.ColumnCount; j++)
                {
                    if(Math.Abs(nodalFreedom[i,j]) > 0.0001)
                    {
                        activeDoF++;
                    }
                }   
            }
            // rank iteration ignored from page 21
            var globalStiffness = Matrix<double>.Build.Dense(activeDoF,activeDoF); // initialize KK matrix
            //setup element stiffness matrix kl (local)\
            var kl0 = BuildElementInLocalStiffness(0);
            var c0 = TransformLocalToGlobal(0);
            var kg0 = BuildElementInGlobalStiffness(kl0,c0);
            var g0 = GetSteeringVector((int)elementConnections[0,0], (int)elementConnections[0, 1]);
            globalStiffness = Form_GSM_KK(kg0, g0, globalStiffness);

            var e2_kg = BuildElementInGlobalStiffness(BuildElementInLocalStiffness(1), TransformLocalToGlobal(1));
            var g2 = GetSteeringVector((int)elementConnections[0, 0], (int)elementConnections[0, 1]);
            globalStiffness = Form_GSM_KK(e2_kg, g2, globalStiffness);

            var e3_kg = BuildElementInGlobalStiffness(BuildElementInLocalStiffness(2), TransformLocalToGlobal(2));
            var g3 = GetSteeringVector((int)elementConnections[0, 0], (int)elementConnections[0, 1]);
            globalStiffness = Form_GSM_KK(e2_kg, g2, globalStiffness);


            Console.WriteLine(globalStiffness.ToMatrixString());
            Console.ReadKey();

        }
        // current GSM is 3x3 - theory says 6x6? investigate further
        //resultant KK matrix is really only valid for element zero
        private static Matrix<double> Form_GSM_KK(Matrix<double> kg, Matrix<double> g, Matrix<double> GSM)
        {
            for (int i = 0; i < activeDoF; i++)
            {
                if (Math.Abs(g[i, 0]) > 0.0001)
                    for (int j = 0; j < activeDoF; j++)
                    {
                        if (Math.Abs(g[j, 0]) > 0.0001)
                        {
                            var gi =i;
                            var gj =j;
                            var val = kg[i, j];
                            var gsmVal = GSM[gi, gj];
                            GSM[gi, gj] = gsmVal+val;
                        }
                            
                    }
            }
            return GSM;
        }


        private static Matrix<double> GetSteeringVector(int i, int j)
        {
            var g = Matrix<double>.Build.Dense(4,1);
            g[0, 0] = nodalFreedom[i, 0];
            g[1, 0] = nodalFreedom[i, 1];
            g[2, 0] = nodalFreedom[j, 0];
            g[3, 0] = nodalFreedom[j, 1];

            return g;
        }

        private static Matrix<double> BuildElementInGlobalStiffness(Matrix<double> kl, Matrix<double> C)
        {
            var kg = C.Multiply(kl).Multiply(C.Transpose()); 
            return kg;
        }

        private static Matrix<double> TransformLocalToGlobal(int i)
        {
            var n1 = elementConnections[i, 0];
            var n2 = elementConnections[i, 1];
            var x1 = nodeCoordinates[(int)n1-1, 0];
            var y1 = nodeCoordinates[(int)n1-1, 1];
            var x2 = nodeCoordinates[(int)n2-1, 0];
            var y2 = nodeCoordinates[(int)n2-1, 1];
            double theta;
            if(x2-x1 == 0)
                if (y2 > y1)
                    theta = 2*Math.Atan(1);
                else
                    theta = -2*Math.Atan(1);
            else
                theta= Math.Atan2(y2 - y1, x2 - x1);
            var C = Matrix<double>.Build.Dense(4, 4); // probably should get rid of magic #

            C[0, 0] = Math.Cos(theta);
            C[0, 1] = -Math.Sin(theta);
            C[0, 2] = 0;
            C[0, 3] = 0;

            C[1, 0] = Math.Sin(theta);
            C[1, 1] = Math.Cos(theta);
            C[1, 2] = 0;
            C[1, 3] = 0;

            C[2, 0] = 0;
            C[2, 1] = 0;
            C[2, 2] = Math.Cos(theta);
            C[2, 3] = -Math.Sin(theta);

            C[3, 0] = 0;
            C[3, 1] = 0;
            C[3, 2] = Math.Sin(theta);
            C[3, 3] = Math.Cos(theta);

            return C;

        }
        
        private static Matrix<double> BuildElementInLocalStiffness(int i)
        {
            var n1 = elementConnections[i, 0];
            var n2 = elementConnections[i, 1];
            var x1 = nodeCoordinates[(int)n1-1, 0];
            var y1 = nodeCoordinates[(int)n1-1, 1];
            var x2 = nodeCoordinates[(int)n2-1, 0];
            var y2 = nodeCoordinates[(int)n2-1, 1];
            var Length = Math.Sqrt(Math.Pow(x2 - x1, 2) + Math.Pow(y2-y1, 2));
            var E = elementProperties[i,0];
            var A = elementProperties[i, 1];
            var kl = Matrix<double>.Build.Dense(4, 4); // probably should get rid of magic #
            var AE_L = A*E/Length;

            kl[0, 0] = AE_L;
            kl[0, 1] = 0;
            kl[0, 2] = -AE_L;
            kl[0, 3] = 0;

            kl[1, 0] = 0;
            kl[1, 1] = 0;
            kl[1, 2] = 0;
            kl[1, 3] = 0;

            kl[2, 0] = -AE_L;
            kl[2, 1] = 0;
            kl[2, 2] = AE_L;
            kl[2, 3] = 0;

            kl[3, 0] = 0;
            kl[3, 1] = 0;
            kl[3, 2] = 0;
            kl[3, 3] = 0;

            return kl;

        }

        private static void InitializeNodeElementMatrices()
        {
            //node XY coordinates
            nodeCoordinates[0, 0] = 0;
            nodeCoordinates[0, 1] = 0;

            nodeCoordinates[1, 0] = 4000;
            nodeCoordinates[1, 1] = 0;

            nodeCoordinates[2, 0] = 4000;
            nodeCoordinates[2, 1] = 6000;

            // element connection matrix
            elementConnections[0, 0] = 1;
            elementConnections[0, 1] = 2;

            elementConnections[1, 0] = 2;
            elementConnections[1, 1] = 3;

            elementConnections[2, 0] = 1;
            elementConnections[2, 1] = 3;
            // element property matrix (E,A)
            elementProperties[0, 0] = 200000;
            elementProperties[0, 1] = 2300;

            elementProperties[1, 0] = 200000;
            elementProperties[1, 1] = 2300;

            elementProperties[2, 0] = 200000;
            elementProperties[2, 1] = 2300;
        }

        private static void InitBCMatrix()
        {
            // 0 means constrained, 1 is free.
            nodalFreedom[0, 0] = 0;
            nodalFreedom[0, 1] = 0;

            nodalFreedom[1, 0] = 1;
            nodalFreedom[1, 1] = 0;

            nodalFreedom[2, 0] = 1;
            nodalFreedom[2, 1] = 1;


            appliedLoad[0, 0] = 0;
            appliedLoad[0, 1] = 0;

            appliedLoad[1, 0] = 0;
            appliedLoad[1, 1] = 0;

            appliedLoad[2, 0] = 1200;
            appliedLoad[2, 1] = 0;

        }

    }
}
