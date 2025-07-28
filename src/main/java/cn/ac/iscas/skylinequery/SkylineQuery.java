package cn.ac.iscas.skylinequery; //skylinequery.SkylineQuery


import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import cn.ac.iscas.utils.RunningTimeCounter;
import cn.ac.iscas.utils.Util;

import static cn.ac.iscas.secretsharing.AdditiveSecretSharing.*;

public class SkylineQuery {

    public static class Vertex {
        public BigInteger[] data;

        public Vertex(int m) {
            data = new BigInteger[m];
        }

        public Vertex(BigInteger[] data) {
            this.data = data;
        }
    }

    public static class Point {
        public BigInteger id;
        public BigInteger[] data;
        public List<Vertex> vertices;

        // public BigInteger label;

        public Point(int m) {
            data = new BigInteger[m];
        }

        public Point(BigInteger id, BigInteger[] data) {
            this.id = id;
            this.data = data;
        }

        public Point(BigInteger id, BigInteger[] data, List<Vertex> vertices) {
            this.id = id;
            this.data = data;
            this.vertices = vertices;
        }
    }

    public static BigInteger[] secureNEuclideanDistance(PartyID partyID, Point[] points, Point queryPoint,
                                                        MultiplicationTriple triple, BigInteger mod, BufferedReader reader, PrintWriter writer) throws IOException {

        int num = points.length;
        int m = queryPoint.data.length;

        BigInteger[] diffis = new BigInteger[num * m];
        for (int i = 0; i < num; i++) {
            int offset = i * m;
            for (int j = 0; j < m; j++) {
                diffis[offset + j] = points[i].data[j].subtract(queryPoint.data[j]).mod(mod);
            }
        }
        BigInteger[] tis = multiplyS(partyID, diffis, diffis, triple, mod, reader, writer);

        BigInteger[] distanceis = new BigInteger[num];
        for (int i = 0; i < num; i++) {
            int offset = i * m;
            distanceis[i] = BigInteger.ZERO;

            for (int j = 0; j < m; j++) {
                distanceis[i] = distanceis[i].add(tis[offset + j]).mod(mod);
            }
        }

        return distanceis;
    }

    public static BigInteger[] secureSum(Point[] points, BigInteger mod) throws IOException {

        int num = points.length;
        int d = points[0].data.length;

        BigInteger[] sumis = new BigInteger[num];
        for (int i = 0; i < num; i++) {
            sumis[i] = BigInteger.ZERO;

            for (int j = 0; j < d; j++) {
                sumis[i] = sumis[i].add(points[i].data[j]).mod(mod);
            }
        }

        return sumis;
    }

    public static Point[] secureBaselineSQ(PartyID partyID, Point[] points, Point[] queryPoints, int q,
                                           MultiplicationTriple triple, RandomNumberTuple rTuple, BigInteger mod,
                                           BufferedReader reader, PrintWriter writer) throws IOException {

        Point[] pointsMapis = new Point[points.length];
        for(int i = 0; i < points.length; i++){
            pointsMapis[i] = new Point(points[0].data.length - 2 + q);
        }

        for(int i = 0; i < q; i++){
            BigInteger[] distanceis = secureNEuclideanDistance(partyID, points, queryPoints[i], triple, mod, reader, writer);
            for(int j = 0; j < points.length; j++){
                pointsMapis[j].data[i] = distanceis[j];
            }
        }

        for(int i = q; i < points[0].data.length - 2 + q; i++){
            for(int j = 0; j < points.length; j++){
                pointsMapis[j].data[i] = points[j].data[i-(q-2)];
            }
        }

        for(int i = 0; i < points.length; i++){
            pointsMapis[i].id = points[i].id;
        }

        BigInteger[] colSumis = secureSum(pointsMapis, mod);

        List<Point> skylinePoints = secureBaselineSQCore(partyID, points, pointsMapis, colSumis, null, triple, rTuple, mod, reader, writer);

        int skylineLen = skylinePoints.size();
        Point[] resulti = new Point[skylineLen];
        for(int i = 0; i < skylineLen; i++){
            resulti[i] = skylinePoints.get(i);
        }

        return resulti;
    }

    private static List<Point> secureBaselineSQCore(PartyID partyID, Point[] points, Point[] pointsMap, BigInteger[] colSum,
                                                    BigInteger[] labels,  MultiplicationTriple triple, RandomNumberTuple rTuple, BigInteger mod,
                                                    BufferedReader reader, PrintWriter writer) throws IOException {

        boolean labelIsNull = (labels == null);
        BigInteger MAX_SUM = mod.divide(BigInteger.TWO).subtract(BigInteger.ONE);

        int num = pointsMap.length;
        int m = pointsMap[0].data.length;
        int d = points[0].data.length;

        //skyline query
        List<Point> skylinePoints = new ArrayList<>();
        boolean isStop = false;
        while (!isStop) {

            int count = 0;
            while (count < 1) {
                int len = num - count;

                while (len > 1) {
                    int offset = (len % 2 == 0) ? count : count + 1;

                    int subLen = len / 2;
                    BigInteger[] leftis = Arrays.copyOfRange(colSum, offset, offset + subLen);
                    BigInteger[] rightis = Arrays.copyOfRange(colSum, offset + subLen, offset + 2 * subLen);

                    BigInteger[] cmpis = secureComparision(partyID, leftis, rightis, triple, rTuple, mod, reader, writer); // <bool(a < b)>

                    int tSize;
                    if (labelIsNull)
                        tSize = (2 + m + d) * subLen; // ids | distances | pointsMap | points
                    else
                        tSize = (2 + m + d) * subLen + subLen; // ids | distances | points | labels

                    BigInteger[] t1is = new BigInteger[tSize]; // <bool(a < b)>
                    BigInteger[] t2is = new BigInteger[tSize]; // <a - b>
                    for (int i = 0; i < subLen; i++) {
                        int lIndex = offset + i, rIndex = lIndex + subLen;

                        // ids
                        t1is[i] = cmpis[i];
                        t2is[i] = points[lIndex].id.subtract(points[rIndex].id).mod(mod);

                        // colSum
                        int dIndex = subLen + i;
                        t1is[dIndex] = cmpis[i];
                        t2is[dIndex] = colSum[lIndex].subtract(colSum[rIndex]).mod(mod);

                        // pointsMap
                        int pmIndex = 2 * subLen + i * m;
                        for (int j = 0; j < m; j++) {
                            t1is[pmIndex + j] = cmpis[i];
                            t2is[pmIndex + j] = pointsMap[lIndex].data[j].subtract(pointsMap[rIndex].data[j]).mod(mod);
                        }

                        // points
                        int pIndex = (2 + m) * subLen + i * d;
                        for (int j = 0; j < d; j++) {
                            t1is[pIndex + j] = cmpis[i];
                            t2is[pIndex + j] = points[lIndex].data[j].subtract(points[rIndex].data[j]).mod(mod);
                        }

                        // label
//                        if (!labelIsNull) {
//                            int labelIndex = (2 + m) * subLen + i;
//                            t1is[labelIndex] = cmpis[i];
//                            t2is[labelIndex] = labels[lIndex].subtract(labels[rIndex]).mod(mod);
//                        }
                    }

                    BigInteger[] mulis = multiplyS(partyID, t1is, t2is, triple, mod, reader, writer); // <bool(a < b)> * <a - b>

                    //<t> = A[left] = <a>, A[left] = <b> + <bool(a < b)> * <a - b>, A[right] = <t> + <b> - A[left]
                    for (int i = 0; i < subLen; i++) {
                        int lIndex = offset + i, rIndex = lIndex + subLen;

                        // ids
                        BigInteger[] ti = conditionSwap(mulis[i], points[lIndex].id, points[rIndex].id, mod);
                        points[lIndex].id = ti[0];
                        points[rIndex].id = ti[1];

                        // colSum
                        int dIndex = subLen + i;
                        ti = conditionSwap(mulis[dIndex], colSum[lIndex], colSum[rIndex], mod);
                        colSum[lIndex] = ti[0];
                        colSum[rIndex] = ti[1];

                        // pointsMap
                        int pmIndex = 2 * subLen + i * m;
                        for (int j = 0; j < m; j++) {
                            ti = conditionSwap(mulis[pmIndex + j], pointsMap[lIndex].data[j], pointsMap[rIndex].data[j], mod);
                            pointsMap[lIndex].data[j] = ti[0];
                            pointsMap[rIndex].data[j] = ti[1];
                        }

                        // points
                        int pIndex = (2 + m) * subLen + i * d;
                        for (int j = 0; j < d; j++) {
                            ti = conditionSwap(mulis[pIndex + j], points[lIndex].data[j], points[rIndex].data[j], mod);
                            points[lIndex].data[j] = ti[0];
                            points[rIndex].data[j] = ti[1];
                        }

                        // label
//                        if (!labelIsNull) {
//                            int labelIndex = (2 + m) * subLen + i;
//                            ti = conditionSwap(mulis[labelIndex], labels[lIndex], labels[rIndex], mod);
//                            labels[lIndex] = ti[0];
//                            labels[rIndex] = ti[1];
//                        }
                    }

                    len = (len % 2 == 0) ? subLen : subLen + 1;
                }
                count++;
            }
            int minIndex = 0;

            //BigInteger[] colSumMin = new BigInteger[]{colSum[minIndex]};
            //BigInteger[] MAXis = new BigInteger[]{shareConstant(partyID, MAX_SUM)};
            BigInteger isStopis = secureComparision(partyID, colSum[minIndex], shareConstant(partyID, MAX_SUM), triple, rTuple, mod, reader, writer);
            //BigInteger[] isStopis = secureComparision(partyID, colSumMin, MAXis, triple, rTuple, mod, reader, writer);
            Util.writeBigInteger(isStopis, writer);
            BigInteger isStopis2 = Util.readBigInteger(reader);
            isStop = (isStopis.add(isStopis2).mod(mod).compareTo(BigInteger.ZERO) == 0);
            if (isStop) {
                break;
            }

            Point skylinePoint = points[minIndex];
            BigInteger[] data = new BigInteger[d];
            for (int i = 0; i < d; i++) {
                data[i] = skylinePoint.data[i];
            }
            skylinePoints.add(new Point(skylinePoint.id,data));
            Point skylinePointMap = pointsMap[minIndex];

            //colSum[minIndex] = shareConstant(partyID, MAX_SUM);

            BigInteger flag = shareConstant(partyID, BigInteger.ZERO);
            BigInteger[] isFirstSky = new BigInteger[num];
            BigInteger[] isEqual = new BigInteger[num];

            BigInteger[] t1i = new BigInteger[m * num]; //
            BigInteger[] t2i = new BigInteger[m * num]; // pointsMap
            for(int i = 0; i < num; i++){
                isEqual[i] = shareConstant(partyID, BigInteger.ONE).subtract(secureComparision(partyID, colSum[minIndex], colSum[i], triple, rTuple, mod, reader, writer));
                isFirstSky[i] = multiply(partyID, isEqual[i], shareConstant(partyID, BigInteger.ONE).subtract(flag).mod(mod), triple, mod, reader, writer);
                flag = flag.add(isFirstSky[i]).mod(mod);

                for (int j = 0; j < m; j++) {
                    t1i[i * m + j] = skylinePointMap.data[j];
                    t2i[i * m + j] = pointsMap[i].data[j];
                }
            }

            BigInteger[] cmpis = secureComparision(partyID, t2i, t1i, triple, rTuple, mod, reader, writer);

            // 1 - (b<a) = (a<=b)
            for(int i = 0; i < num * m; i++) {
                cmpis[i] = shareConstant(partyID, BigInteger.ONE).subtract(cmpis[i]).mod(mod);
            }
            BigInteger[] isDomis = new BigInteger[num];
            for (int i = 0; i < num; i++) {
                BigInteger isDomi = shareConstant(partyID, BigInteger.ONE);
                for (int j = 0; j < m; j++){
                    isDomi = multiply(partyID, isDomi, cmpis[i * m + j], triple, mod, reader, writer);
                }
                isDomi = multiply(partyID, isDomi, shareConstant(partyID, BigInteger.ONE).subtract(isEqual[i]).mod(mod), triple, mod, reader, writer);
                isDomis[i] = isDomi.add(isFirstSky[i]).mod(mod);
            }

            t1i = new BigInteger[num * 2]; // (a<=b) | 1-(a<=b)
            t2i = new BigInteger[num * 2]; //   MAX  | colSum
            for (int i = 0; i < num; i++) {
                t1i[i] = isDomis[i];
                t2i[i] = shareConstant(partyID, MAX_SUM);
            }
            for (int i = 0; i < num; i++) {
                t1i[num + i] = shareConstant(partyID, BigInteger.ONE).subtract(isDomis[i]).mod(mod);
                t2i[num + i] = colSum[i];
            }
            //(a<=b)*MAX | (1-(a<=b))*colSum
            BigInteger[] mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);
            //(a<=b)*MAX + (1-(a<=b))*colSum
            for(int i = 0; i < num; i++) {
                colSum[i] = mulis[i].add(mulis[num + i]).mod(mod);
            }

        }
        return skylinePoints;
    }

    // SQ1
    public static Point[] secureSQ1(PartyID partyID, VC[] vcs, Point[] queryPoints, int q, int d,
                                           MultiplicationTriple triple, RandomNumberTuple rTuple, BigInteger mod,
                                           BufferedReader reader, PrintWriter writer) throws IOException {

        int vcNum = vcs.length;
        int vcSize = vcs[0].points.length;
        int m = queryPoints[0].data.length;

        //  bool( low_i <= q_i < high_i ) = ( 1 - bool( q_i <= low_i ) * bool( q_i < high_i)
        BigInteger[] t1i = new BigInteger[vcNum * m * 2]; //   q_1, ..., q_m   | q_1, ..., q_m
        BigInteger[] t2i = new BigInteger[vcNum * m * 2]; // low_1, ..., low_m | high_1, ..., high_m
        for (int i = 0; i < vcNum; i++) {
            int index = i * m * 2;
            for (int j = 0; j < m; j++) {
                t1i[index + j] = queryPoints[0].data[j];
                t2i[index + j] = vcs[i].low.data[j];

                t1i[index + m + j] = queryPoints[0].data[j];
                t2i[index + m + j] = vcs[i].high.data[j];
            }
        }

        BigInteger[] cmpis = secureComparision(partyID, t1i, t2i, triple, rTuple, mod, reader, writer);

        // bool( low_i <= q_i ) = 1 - bool( q_i < low_i )
        for (int i = 0; i < vcNum; i++) {
            int index = i * m * 2;

            for (int j = 0; j < m; j++) {
                cmpis[index + j] = shareConstant(partyID, BigInteger.ONE).subtract(cmpis[index + j]).mod(mod);
            }
        }

        // PROD( bool(low_i <= q_i) * bool(q_i < high_i) )ji
        t1i = new BigInteger[vcNum * m]; // bool(low_i <= q_i)
        t2i = new BigInteger[vcNum * m]; // bool(q_i < high_i)
        for (int i = 0; i < vcNum; i++) {
            int index1 = i * m;
            int index2 = i * m * 2;

            for (int j = 0; j < m; j++) {
                t1i[index1 + j] = cmpis[index2 + j];
                t2i[index1 + j] = cmpis[index2 + m + j];
            }
        }
        BigInteger[] mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        t1i = new BigInteger[vcNum];
        t2i = new BigInteger[vcNum];
        for (int i = 0; i < vcNum; i++) {
            int index = i * 2;

            t1i[i] = mulis[index];
            t2i[i] = mulis[index + 1];
        }

        BigInteger[] alphais = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        SkylineQuery.Point[][] pDatas = new SkylineQuery.Point[vcNum][];
        for (int i = 0; i < vcNum; i++) {
            pDatas[i] = vcs[i].points;
        }
        int vNum = vcs[0].points[0].vertices.size();

        SkylineQuery.Point[] pointis = new SkylineQuery.Point[vcSize]; // 候选点集

        obliGetData(partyID, pointis, vcNum, vcSize, d, vNum, alphais, pDatas,
                triple, mod, reader, writer);

        Point[] resulti = secureSQ1Core(partyID, pointis, queryPoints, triple, rTuple, mod, reader, writer);

        return resulti;
    }

    private static Point[] secureSQ1Core(PartyID partyID, Point[] points, Point[] queryPoints, MultiplicationTriple triple, RandomNumberTuple rTuple, BigInteger mod,
                                                    BufferedReader reader, PrintWriter writer) throws IOException {
        int pNum = points.length;
        int d = points[0].data.length;
        int v = points[0].vertices.size();

        BigInteger[] t1i = new BigInteger[pNum * v * 2]; //   L1x   L1y  L2x L2y ...  Lvx Lvy  | L1x L1y L2x L2y ...Lvx Lvy | ...
        BigInteger[] t2i = new BigInteger[pNum * v * 2]; //   qv1y qv1x qv2y qv2x ...qvvy qvvx | ...
        for (int i = 0; i < pNum; i++) {
            int index = i * v * 2;
            for (int j = 0; j < v; j++) {
                t1i[index + j * 2] = points[i].vertices.get((j + 1) % v).data[0].subtract(points[i].vertices.get(j).data[0]).mod(mod);
                t1i[index + j * 2 + 1] = points[i].vertices.get((j + 1) % v).data[1].subtract(points[i].vertices.get(j).data[1]).mod(mod);

                t2i[index + j * 2] = queryPoints[0].data[1].subtract(points[i].vertices.get(j).data[1]).mod(mod);
                t2i[index + j * 2 + 1] = queryPoints[0].data[0].subtract(points[i].vertices.get(j).data[0]).mod(mod);
            }
        }
        BigInteger[] mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        BigInteger[] subis = new BigInteger[pNum * v];

//        BigInteger newMod = BigInteger.probablePrime(((int) Util.log2(mod.doubleValue()))*2, new Random());
        for (int i = 0; i < pNum; i++) {
            int index = i * v;
            for (int j = 0; j < v; j++) {
//                subis[index + j] = mulis[index * 2 + j * 2].subtract(mulis[index * 2 + j * 2 + 1]).mod(mod); // L1x * qv1y - L1y * qv1x
                subis[index + j] = mulis[index * 2 + j * 2 + 1].subtract(mulis[index * 2 + j * 2]).mod(mod); // L1y * qv1x - L1x * qv1y
//                subis[index + j] = mulis[index * 2 + j * 2 + 1].subtract(mulis[index * 2 + j * 2]).mod(newMod); // L1y * qv1x - L1x * qv1y
            }
        }

        BigInteger[] cmpis = secureComparisionSub1(partyID, subis, triple, rTuple, mod, reader, writer);

        BigInteger[] cmpis1 = new BigInteger[pNum];
        BigInteger[][] cmpisP = new BigInteger[v][pNum];
        Arrays.fill(cmpis1, shareConstant(partyID, BigInteger.ONE));
        for (int i = 0; i < v; i++) {
            for (int j = 0; j < pNum; j++) {
                cmpisP[i][j] = cmpis[j * v + i];
            }
        }
        for (int i = 0; i < v; i++) {
            cmpis1 = multiplyS(partyID, cmpis1, cmpisP[i], triple, mod, reader, writer);
        }

        BigInteger[]  t4i = new BigInteger[cmpis.length]; //   0 0 0 0 | ...
        for (int i = 0; i < cmpis.length; i++){
            t4i[i] = shareConstant(partyID, BigInteger.ZERO);
        }

        BigInteger[] cmpis2 = secureEqual(partyID, subis, t4i, triple, rTuple, mod, reader, writer);

        BigInteger[] cmpis3 = new BigInteger[pNum];
        Arrays.fill(cmpis3, shareConstant(partyID, BigInteger.ONE));
        for (int i = 0; i < v; i++) {
            for (int j = 0; j < pNum; j++) {
                cmpisP[i][j] = cmpis2[j * v + i];
            }
        }
        for (int i = 0; i < v; i++) {
            cmpis3 = multiplyS(partyID, cmpis3, cmpisP[i], triple, mod, reader, writer);
        }
        for (int i = 0; i < cmpis3.length; i++){
            cmpis3[i] = shareConstant(partyID, BigInteger.ONE).subtract(cmpis3[i]).mod(mod);
        }

        BigInteger[] resulti = multiplyS(partyID, cmpis1, cmpis3, triple, mod, reader, writer);

        t1i = new BigInteger[pNum * (1 + d)]; //   id   | pointData
        t2i = new BigInteger[pNum * (1 + d)]; // result |  result
        for (int i = 0; i < pNum; i++) {
            int index1 = i * (1 + d);

            t1i[index1] = points[i].id;
            System.arraycopy(points[i].data, 0, t1i, index1 + 1, d);

            Arrays.fill(t2i, index1, index1 + 1 + d, resulti[i]);
        }
        mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        Point[] skylineListi = new Point[pNum];
        BigInteger id;
        BigInteger[] data;
        for (int i = 0; i < pNum; i++) {
            int index1 = i * (1 + d);
            id = mulis[index1];
            data = new BigInteger[d];
            for (int l = 0; l < d; l++) {
                data[l] = mulis[index1 + 1 + l];
            }
            skylineListi[i] = new Point(id, data);
        }

        return skylineListi;
    }

    private static Point[] secureSQ1Core2(PartyID partyID, Point[] points, Point[] queryPoints, MultiplicationTriple triple, RandomNumberTuple rTuple, BigInteger mod,
                                         BufferedReader reader, PrintWriter writer) throws IOException {
        int pNum = points.length;
        int d = points[0].data.length;
        int v = points[0].vertices.size();

        BigInteger[] t0i = new BigInteger[pNum * v * 2]; //   v1x  v1y  v2x v2y ...  vvx vvy  | v1x v1y v2x v2y ...vvx vvy | ...
        for (int i = 0; i < pNum; i++) {
            int index = i * v * 2;
            for (int j = 0; j < v; j++) {
                t0i[index + j * 2] = points[i].vertices.get(j).data[0];
                t0i[index + j * 2 + 1] = points[i].vertices.get(j).data[1];
            }
        }
        BigInteger[] t2i = new BigInteger[pNum * v]; //   v1y  v2y ...  vvy | v1y  v2y ...  vvy | ...
        for (int i = 0; i < pNum; i++) {
            int index = i * v;
            for (int j = 0; j < v; j++) {
                t2i[index + j] = t0i[index * 2 + j * 2 + 1];
            }
        }

        BigInteger[] t1i = new BigInteger[pNum * v * 2]; //   L1x   L1y  L2x L2y ...  Lvx Lvy  | L1x L1y L2x L2y ...Lvx Lvy | ...
//        BigInteger[] t2i = new BigInteger[pNum * v * 2]; //   qv1y qv1x qv2y qv2x ...qvvy qvvx | ...
        for (int i = 0; i < pNum; i++) {
            int index = i * v * 2;
            for (int j = 0; j < v; j++) {
                t1i[index + j * 2] = points[i].vertices.get((j + 1) % v).data[0].subtract(points[i].vertices.get(j).data[0]).mod(mod);
                t1i[index + j * 2 + 1] = points[i].vertices.get((j + 1) % v).data[1].subtract(points[i].vertices.get(j).data[1]).mod(mod);

//                t2i[index + j * 2] = queryPoints[0].data[1].subtract(points[i].vertices.get(j).data[1]).mod(mod);
//                t2i[index + j * 2 + 1] = queryPoints[0].data[0].subtract(points[i].vertices.get(j).data[0]).mod(mod);
            }
        }

//        BigInteger[] t3i = new BigInteger[pNum * v]; //   v1y  v2y ...  vvy | v1y  v2y ...  vvy | ...
//        System.arraycopy(t2i, 0, t3i, 0, pNum * v);

        BigInteger[] t4i = new BigInteger[pNum * v]; //   v2y v3y v4y v1y ...
        for (int i = 0; i < pNum; i++) {
            int index = i * v;
            for (int j = 0; j < v; j++) {
                t4i[index + j] = t2i[index + ((j + 1) % v)];
            }
        }
        BigInteger[] equis0 = secureComparision(partyID, t2i, t4i, triple, rTuple, mod, reader, writer);
        BigInteger[] equis1 = secureComparision(partyID, t4i, t2i, triple, rTuple, mod, reader, writer);
        BigInteger[] equis = new BigInteger[pNum * v];
        for (int i = 0; i < pNum * v; i++) {
            equis[i] = equis0[i].add(equis1[i]).mod(mod);
        }

        BigInteger[] t5i = new BigInteger[pNum * v * 2]; //   v1y v2y v3y v4y ...  | v2y v3y v4y v1y ...

        System.arraycopy(t2i, 0, t5i, 0, pNum * v);
        for (int i = 0; i < pNum; i++) {
            int index = i * v;
            for (int j = 0; j < v; j++) {
                t5i[pNum * v + index + j] = t2i[index + ((j + 1) % v)];
            }
        }

        int subLen = pNum * v;
        BigInteger[] leftis = Arrays.copyOfRange(t5i, 0, subLen);
        BigInteger[] rightis = Arrays.copyOfRange(t5i, subLen, 2 * subLen);

        BigInteger[] mmis = secureComparision(partyID, leftis, rightis, triple, rTuple, mod, reader, writer); // <bool(a < b)>

        BigInteger[] t6is = new BigInteger[subLen]; // <bool(a < b)>
        BigInteger[] t7is = new BigInteger[subLen]; // <a - b>
        for (int i = 0; i < subLen; i++) {
            int rIndex = i + subLen;

            t6is[i] = mmis[i];
            t7is[i] = t5i[i].subtract(t5i[rIndex]).mod(mod);
        }

        BigInteger[] mulis = multiplyS(partyID, t6is, t7is, triple, mod, reader, writer); // <bool(a < b)> * <a - b>

        // <t> = A[left] = <a>, A[left] = <b> + <bool(a < b)> * <a - b>, A[right] = <t> + <b> - A[left]
        for (int i = 0; i < subLen; i++) {
            int rIndex = i + subLen;

            BigInteger[] ti = conditionSwap(mulis[i], t5i[i], t5i[rIndex], mod);
            t5i[i] = ti[0];
            t5i[rIndex] = ti[1];
        }

        //querypoint.y与min(y1,y2)/max(y1,y2)
        BigInteger[] qyi = new BigInteger[subLen]; // qy
        Arrays.fill(qyi, queryPoints[0].data[1]);

        BigInteger[] Lmincmpis = secureComparision(partyID, qyi, Arrays.copyOfRange(t5i, 0, subLen), triple, rTuple, mod, reader, writer);
        for (int i = 0; i < subLen; i++) {
            Lmincmpis[i] = shareConstant(partyID, BigInteger.ONE).subtract(Lmincmpis[i]).mod(mod);
        }
        BigInteger[] Lmaxcmpis = secureComparision(partyID, qyi, Arrays.copyOfRange(t5i, subLen, subLen * 2), triple, rTuple, mod, reader, writer);

        // x=(qy-v1y)*(v2x-v1x)/(v2y-v1y) + v1x = (qy-v1y)*(v2x-v1x)*(v2y-v1y)^-1 + v1x
        BigInteger[] t9i = new BigInteger[subLen];
        for (int i = 0; i < subLen; i++) {
            t9i[i] = qyi[i].subtract(t2i[i]).mod(mod);
        }
        BigInteger[] txxi = new BigInteger[subLen];
        for (int i = 0; i < subLen; i++){
            txxi[i] = t1i[i * 2];
        }
        BigInteger[] tyyi = new BigInteger[subLen];
        for (int i = 0; i < subLen; i++){
            tyyi[i] = t1i[i * 2 + 1];
        }
        for (int i = 0; i < subLen; i++){
            int power = mod.intValue()-2;
            BigInteger base = tyyi[i];
            tyyi[i] = shareConstant(partyID, BigInteger.ONE);
            while (power > 0){
                if (power % 2 == 1){
                    tyyi[i] = multiply(partyID, tyyi[i], base, triple, mod, reader, writer);
                }
                base = multiply(partyID, base, base, triple, mod, reader, writer);
                power = power >> 1;
            }
        }
        BigInteger[] yxmulis = multiplyS(partyID, t9i, txxi, triple, mod, reader, writer);
        BigInteger[] yxymulis = multiplyS(partyID, yxmulis, tyyi, triple, mod, reader, writer);
        BigInteger[] xi = new BigInteger[subLen];
        for (int i = 0; i < subLen; i++){
            xi[i] = yxymulis[i].add(t0i[i * 2]).mod(mod);
        }
        BigInteger[] qxi = new BigInteger[subLen]; // qx
        Arrays.fill(qxi, queryPoints[0].data[0]);
        BigInteger[] qxxcmpis = secureComparision(partyID, qxi, xi, triple, rTuple, mod, reader, writer);

        qxxcmpis = multiplyS(partyID, qxxcmpis, Lmincmpis, triple, mod, reader, writer);
        qxxcmpis = multiplyS(partyID, qxxcmpis, Lmaxcmpis, triple, mod, reader, writer);
        qxxcmpis = multiplyS(partyID, qxxcmpis, equis, triple, mod, reader, writer);

        BigInteger[] resulti = new BigInteger[pNum];
        for (int i = 0; i < pNum; i++) {
            int index = i * v;
            resulti[i] = shareConstant(partyID, BigInteger.ZERO);
            for (int j = 0; j < v; j++) {
                resulti[i] = resulti[i].add(qxxcmpis[index + j]).mod(mod);
            }
        }

        t1i = new BigInteger[pNum * (1 + d)]; //   id   | pointData
        t2i = new BigInteger[pNum * (1 + d)]; // result |  result
        for (int i = 0; i < pNum; i++) {
            int index1 = i * (1 + d);

            t1i[index1] = points[i].id;
            System.arraycopy(points[i].data, 0, t1i, index1 + 1, d);

            Arrays.fill(t2i, index1, index1 + 1 + d, resulti[i]);
        }
        mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        Point[] skylineListi = new Point[pNum];
        BigInteger id;
        BigInteger[] data;
        for (int i = 0; i < pNum; i++) {
            int index1 = i * (1 + d);
            id = mulis[index1];
            data = new BigInteger[d];
            for (int l = 0; l < d; l++) {
                data[l] = mulis[index1 + 1 + l];
            }
            skylineListi[i] = new Point(id, data);
        }

        return skylineListi;
    }

    private static BigInteger[] conditionSwap(BigInteger muli, BigInteger ai, BigInteger bi, BigInteger mod) {
        BigInteger[] resulti = new BigInteger[2];
        resulti[0] = bi.add(muli).mod(mod);
        resulti[1] = ai.add(bi).subtract(resulti[0]).mod(mod);

        return resulti;
    }

    public static class AG {

        public BigInteger label;
        public Point[] points;
        public BigInteger[] subLabels;

        public AG(int num, int m) {
            points = new Point[num];
            for (int i = 0; i < num; i++) {
                points[i] = new Point(m);
            }

            subLabels = new BigInteger[num];
        }

        public AG(BigInteger label, Point[] points, BigInteger[] subLabels) {
            this.label = label;
            this.points = points;
            this.subLabels = subLabels;
        }
    }

    public static class VG {

        public Point low, high; // 边界
        public Point[] points; // 点
        public BigInteger[] subLabels; // 各点对应的label

        public VG(int num, int m) {
            points = new Point[num];
            for (int i = 0; i < num; i++) {
                points[i] = new Point(m);
            }

            subLabels = new BigInteger[num];
        }

        public VG(Point low, Point high, Point[] points, BigInteger[] subLabels) {
            this.low = low;
            this.high = high;
            this.points = points;
            this.subLabels = subLabels;
        }

    }

    public static class VC {

        public Point low, high; // 边界
        public Point[] points; // 点

        public VC(int num, int m) {
            points = new Point[num];
            for (int i = 0; i < num; i++) {
                points[i] = new Point(m);
            }
        }

        public VC(Point low, Point high, Point[] points) {
            this.low = low;
            this.high = high;
            this.points = points;
        }

    }

    private static void getSelectedData(PartyID partyID, Point[] pointis, BigInteger[] labelis, int num, int size,
                                        int m, BigInteger[] alphais, Point[][] pDatas, BigInteger[][] lDatas,
                                        MultiplicationTriple triple, BigInteger mod, BufferedReader reader, PrintWriter writer) throws IOException {

        int tSize = size * (2 + m);
        BigInteger[] t1i = new BigInteger[num * tSize]; // label |  id   | pointData
        BigInteger[] t2i = new BigInteger[num * tSize]; // alpha | alpha | alpha
        for (int i = 0; i < num; i++) {
            int index1 = i * tSize;

            for (int j = 0; j < size; j++) {
                int index2 = index1 + j * (2 + m);

                t1i[index2] = lDatas[i][j];
                t1i[index2 + 1] = pDatas[i][j].id;
                System.arraycopy(pDatas[i][j].data, 0, t1i, index2 + 2, m);
            }

            Arrays.fill(t2i, index1, index1 + tSize, alphais[i]);
        }
        BigInteger[] mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        for (int i = 0; i < size; i++) {
            labelis[i] = BigInteger.ZERO;

            pointis[i] = new Point(m);
            pointis[i].id = BigInteger.ZERO;
            pointis[i].data = new BigInteger[m];
            for (int j = 0; j < m; j++) {
                pointis[i].data[j] = BigInteger.ZERO;
            }
        }

        for (int i = 0; i < num; i++) {
            int index1 = i * tSize;

            for (int j = 0; j < size; j++) {
                int index2 = index1 + j * (2 + m);

                labelis[j] = labelis[j].add(mulis[index2]).mod(mod);
                pointis[j].id = pointis[j].id.add(mulis[index2 + 1]).mod(mod);

                for (int l = 0; l < m; l++) {
                    pointis[j].data[l] = pointis[j].data[l].add(mulis[index2 + 2 + l]).mod(mod);
                }
            }
        }
    }

    private static void obliGetData(PartyID partyID, Point[] pointis, int num, int size,
                                        int d, int v, BigInteger[] alphais, Point[][] pDatas,
                                        MultiplicationTriple triple, BigInteger mod, BufferedReader reader, PrintWriter writer) throws IOException {

        int tSize = size * (1 + d + v * 2);
        BigInteger[] t1i = new BigInteger[num * tSize]; //  id   | pointData | vertex
        BigInteger[] t2i = new BigInteger[num * tSize]; // alpha |   alpha   | alpha
        for (int i = 0; i < num; i++) {
            int index1 = i * tSize;

            for (int j = 0; j < size; j++) {
                int index2 = index1 + j * (1 + d + v * 2);

                t1i[index2] = pDatas[i][j].id;
                System.arraycopy(pDatas[i][j].data, 0, t1i, index2 + 1, d);
                for (int l = 0; l < v; l++) {
                    System.arraycopy(pDatas[i][j].vertices.get(l).data, 0, t1i, index2 + 1 + d + 2 * l, 2);
                }
            }

            Arrays.fill(t2i, index1, index1 + tSize, alphais[i]);
        }
        BigInteger[] mulis = multiplyS(partyID, t1i, t2i, triple, mod, reader, writer);

        for (int i = 0; i < size; i++) {

            pointis[i] = new Point(d);
            pointis[i].id = shareConstant(partyID, BigInteger.ZERO);
            pointis[i].data = new BigInteger[d];
            for (int j = 0; j < d; j++) {
                pointis[i].data[j] = shareConstant(partyID, BigInteger.ZERO);
            }
            pointis[i].vertices = new ArrayList<>();
            for (int l = 0; l < v; l++) {
                pointis[i].vertices.add(new Vertex(new BigInteger[]{shareConstant(partyID, BigInteger.ZERO), shareConstant(partyID, BigInteger.ZERO)}));
            }

        }

        for (int i = 0; i < num; i++) {
            int index1 = i * tSize;

            for (int j = 0; j < size; j++) {
                int index2 = index1 + j * (1 + d + v * 2);

                pointis[j].id = pointis[j].id.add(mulis[index2]).mod(mod);

                for (int l = 0; l < d; l++) {
                    pointis[j].data[l] = pointis[j].data[l].add(mulis[index2 + 1 + l]).mod(mod);
                }
                for (int l = 0; l < v; l++) {
                    pointis[j].vertices.get(l).data[0] = pointis[j].vertices.get(l).data[0].add(mulis[index2 + 1 + d + 2 * l]).mod(mod);
                    pointis[j].vertices.get(l).data[1] = pointis[j].vertices.get(l).data[1].add(mulis[index2 + 1 + d + 2 * l + 1]).mod(mod);
                }
            }
        }
    }

    public static void main(String[] args) {

        List<BigInteger> list = new ArrayList<>();
        list.add(BigInteger.valueOf(0L));
        list.add(BigInteger.valueOf(1L));
        list.add(BigInteger.valueOf(2L));

        BigInteger[] array = list.toArray(new BigInteger[] {});
        System.out.println(array);
    }
}

