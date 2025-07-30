package cn.ac.iscas;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.IntStream;

import cn.ac.iscas.kdtree.KDTreeNode;
import cn.ac.iscas.kdtree.KDTreePoint;
import cn.ac.iscas.voronoi.VoronoiWithTinfour;
import com.alibaba.fastjson2.JSON;

import static cn.ac.iscas.secretsharing.AdditiveSecretSharing.*;

import cn.ac.iscas.secretsharing.AdditiveSecretSharing;
import cn.ac.iscas.secretsharing.AdditiveSecretSharing.MultiplicationTriple;
import cn.ac.iscas.secretsharing.AdditiveSecretSharing.RandomNumberTuple;
import cn.ac.iscas.skylinequery.SkylineQuery;
import cn.ac.iscas.skylinequery.SkylineQuery.AG;
import cn.ac.iscas.skylinequery.SkylineQuery.Point;
import cn.ac.iscas.skylinequery.SkylineQuery.Vertex;
import cn.ac.iscas.skylinequery.SkylineQuery.VG;
import cn.ac.iscas.skylinequery.SkylineQuery.VC;
import cn.ac.iscas.utils.RunningTimeCounter;
import cn.ac.iscas.utils.Util;
import com.alibaba.fastjson2.JSONWriter;
import org.tinfour.demo.viewer.backplane.LidarPointSelection;

import static cn.ac.iscas.utils.DataProcessor.*;
import static cn.ac.iscas.utils.QuickSelect.findKSmallestIndexes;

public class TestSSQV1 {

    /**

     * @param args
     * @throws IOException
     */
    public static void user(String[] args) throws IOException {

        int index = 1;
        String ipC1 = args[index++];
        int portC1 = Integer.parseInt(args[index++]);
        String ipC2 = args[index++];
        int portC2 = Integer.parseInt(args[index++]);

        int testType = Integer.parseInt(args[index++]); //testType
        int epoch = Integer.parseInt(args[index++]); //epoch
        String randomSeed = args[index++];
        int n = Integer.parseInt(args[index++]); //dataNumber
        int dataLength = Integer.parseInt(args[index++]); //dataLength
        int d = Integer.parseInt(args[index++]); // dimension
        int q = Integer.parseInt(args[index++]); //queryNumber
        int dt = Integer.parseInt(args[index++]); //dataset

        int l = dataLength * 2 + (int) Util.log2(d) + 2;

        Random random = randomSeed.equals("null") ? new Random() : new Random(Long.parseLong(randomSeed));
        BigInteger mod = BigInteger.probablePrime(l, random);
        MultiplicationTriple[] triples = generateMultiplicationTriples(mod);
        RandomNumberTuple[] tuples = generateRandomNumberTuples(l, mod);
        System.out.println("mod = " + mod);

        BigInteger[][] dataset = null;
        switch (dt){
            case 0: dataset = generateDataset_CORR(n, d, dataLength, random);
                break;
            case 1: dataset = generateDataset_INDE(n, d, dataLength, random);
                break;
            case 2: dataset = generateDataset_ANTI(n, d, dataLength, random);
                break;
            case 3: dataset = getDataset_HOTE(n, d, dataLength, random);
                break;
        }

        BigInteger[][] queryPoints = generateDataset(q, 2, dataLength, random);
        Point[][] pointsSecretShares = null;

        //AG[] ags = null;
        //VG[] vgs = null;
        VC[] vcs = null;
        //AG[][] agsSecrets = null;
        //VG[][] vgsSecrets = null;
        VC[][] vcsSecrets = null;
        //skyline pre-result
        //List<Point> skylinePreResult = new ArrayList<>();
        if (testType == 0) {
            pointsSecretShares = new Point[2][n];
            for (int i = 0; i < n; i++) {
                pointsSecretShares[0][i] = new Point(d);
                pointsSecretShares[1][i] = new Point(d);

                BigInteger[] idSecretShares = randomSplit(dataset[i][d], mod);
                pointsSecretShares[0][i].id = idSecretShares[0];
                pointsSecretShares[1][i].id = idSecretShares[1];

                for (int j = 0; j < d; j++) {
                    BigInteger[] pointSecretShares = randomSplit(dataset[i][j], mod);

                    pointsSecretShares[0][i].data[j] = pointSecretShares[0];
                    pointsSecretShares[1][i].data[j] = pointSecretShares[1];
                }
            }
        } else if (testType == 1) {

            List<Point> pointList= parseDataset(dataset, d);
            List<List<Vertex>> scope = new ArrayList<>();
           // long timeOfIndexConsStart = System.currentTimeMillis();

//            final ExecutorService executor = new ThreadPoolExecutor(16, 17, 60L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<>());
//            final ReentrantLock lock = new ReentrantLock();
//            CountDownLatch latch = new CountDownLatch(pointList.size()); //pointList.size() - 200000
            for (int i = 0; i < pointList.size(); i++) {

                int finalI = i;
//                executor.execute(new Runnable() {
//                    @Override
//                    public void run() {
//                        try{
                List<Point> pointListCopy = new ArrayList<>(pointList);
                List<Point> dominatorsList = getNonspatialDominators2(pointListCopy, pointList.get(finalI), finalI, d);
                // if (dominatorsList.isEmpty()) {
                //     skylinePreResult.add(pointList.get(finalI));
                // }else {
                List<Vertex> verticesList = getVertices(pointList.get(finalI), dominatorsList, finalI, d, dataLength, n, dt);
                scope.add(verticesList);
                //}
//                        } finally {
//                            latch.countDown();
//                        }
//
//                    }
//                });
//                scope.add(verticesList);
            }
//            try {
//                latch.await();
//            } catch (InterruptedException e) {
//                e.printStackTrace();
//            }
//            executor.shutdown();

//            long timeOfIndexCons = System.currentTimeMillis() - timeOfIndexConsStart;
//            System.out.println("Index construction time: " + timeOfIndexCons + " ms");

            //int minNum = scope.get(0).size();
            int maxNum = scope.get(0).size();
            for (int i = 1; i < scope.size(); i++){
//                if(scope.get(i).size() < minNum){
//                    minNum = scope.get(i).size();
//                }
                if(scope.get(i).size() > maxNum){
                    maxNum = scope.get(i).size();
                }
            }
//            int maxNum = 0;
//            for (Map.Entry<Integer, List<Vertex> > entry : scope.entrySet()){
//                if (entry.getValue().size() > maxNum){
//                    maxNum = entry.getValue().size();
//                }
//            }

            for (int i = 0; i < scope.size(); i++){
//                int padNum = random.nextInt(100) % (maxNum - scope.get(i).size() + 1);
                int padNum = maxNum - scope.get(i).size();
                for (int j = 0; j < padNum; j++){
                    scope.get(i).add(scope.get(i).get(scope.get(i).size()-1));
                }
            }

            int minVCSize = Integer.parseInt(args[index++]);
//            int vcNum = (int) Math.ceil(n / minVCSize);

            KDTreePoint[] kdPoints = new KDTreePoint[scope.size()];
            // List<KDTreePoint> kdPoints = new ArrayList<>();

             for (int i = 0; i < n; i++) {
                 if (scope.get(i) == null){
                     continue;
                 }
                 BigInteger id = dataset[i][d];
                 BigInteger[] data = new BigInteger[2];
                 BigInteger maxX = scope.get(i).get(0).data[0];
                 BigInteger maxY = scope.get(i).get(0).data[1];
                 for (int j = 1; j < scope.get(i).size(); j++) {
                     if (scope.get(i).get(j).data[0].compareTo(maxX) > 0) {
                         maxX = scope.get(i).get(j).data[0];
                     }
                     if (scope.get(i).get(j).data[1].compareTo(maxY) > 0) {
                         maxY = scope.get(i).get(j).data[1];
                     }
                 }
                 data[0] = maxX;
                 data[1] = maxY;
                 kdPoints[i] = new KDTreePoint(id, data);
                 //kdPoints.add(new KDTreePoint(id, data));
             }

//            for (int i = 0; i < n; i++) {
//                BigInteger id = dataset[i][d];
//                BigInteger[] data = new BigInteger[2];
//                System.arraycopy(dataset[i], 0, data, 0, 2);
//                kdPoints[i] = new KDTreePoint(id, data);
//            }

            BigInteger[] lb = new BigInteger[] { BigInteger.valueOf(0), BigInteger.valueOf(0) };
            BigInteger[] ub = new BigInteger[] { BigInteger.valueOf(2).pow(dataLength), BigInteger.valueOf(2).pow(dataLength) };
            KDTreeNode kdTree = new KDTreeNode(minVCSize, 2, kdPoints, null, lb, ub);
            List<KDTreeNode> leafNodes = KDTreeNode.getAllLeafPoints(kdTree);

            int vcNum = leafNodes.size();

            List<List<Point>> VCPointLists = new ArrayList<>();
            for (int i = 0; i < leafNodes.size(); i++) {
                List<Point> VCPointList = new ArrayList<>();
                for (int j = 0; j < leafNodes.get(i).points.length; j++) {
                    BigInteger id = leafNodes.get(i).points[j].id;
//                    BigInteger[] data = new BigInteger[d];
//                    List<Vertex> vertices = new ArrayList<>();
//                    for (int k = 0; k < pointList.size(); k++) {
//                        if (id.equals(pointList.get(k).id)) {
//                            data = pointList.get(k).data;
//                            vertices = scope.get(k);
//                            break;
//                        }
//                    }
                    BigInteger[] data = pointList.get(id.intValue()-1).data;
                    List<Vertex> vertices = scope.get(id.intValue()-1);
                    VCPointList.add(new Point(id, data, vertices));
                }
                VCPointLists.add(VCPointList);
            }

            List<List<Point>> VCPointListsTmp0 = new ArrayList<>();
            for (int i = 0; i < VCPointLists.size(); i++) {
                List<Point> pointListTmp = new ArrayList<>(VCPointLists.get(i));
                VCPointListsTmp0.add(pointListTmp);
            }
            for (int i = 0; i < VCPointLists.size(); i++) {
                Point low, high;
                low = new Point(leafNodes.get(i).lb.id, leafNodes.get(i).lb.data);
                high = new Point(leafNodes.get(i).ub.id, leafNodes.get(i).ub.data);
                List<List<Point>> VCPointListsTmp = new ArrayList<>(VCPointListsTmp0);
                VCPointListsTmp.remove(i);
                for (int j = 0; j < VCPointListsTmp.size(); j++) {
                    for (int k = 0; k < VCPointListsTmp.get(j).size(); k++) {
                        if (isIntersect(VCPointListsTmp.get(j).get(k).vertices, low, high)){
                            VCPointLists.get(i).add(VCPointListsTmp.get(j).get(k));
                        }
                    }
                }
            }

            int maxVCSize = VCPointLists.get(0).size();
            for (int i = 1; i < VCPointLists.size(); i++){
                if(VCPointLists.get(i).size() > maxVCSize){
                    maxVCSize = VCPointLists.get(i).size();
                }
            }
            int idAccumulator = n + 1;
            for (int i = 0; i < VCPointLists.size(); i++){
                int padNum = maxVCSize - VCPointLists.get(i).size();

                for (int j = 0; j < padNum; j++){
                    BigInteger dummyId = BigInteger.valueOf(j).add(BigInteger.valueOf(idAccumulator));
                    BigInteger[] dummyData = new BigInteger[d];
                    for (int k = 0; k < d; k++){
                        dummyData[k] = BigInteger.ZERO;
                    }
                    List<Vertex> dummyVertices = new ArrayList<>(maxNum);
                    int tempY = leafNodes.get(i).ub.data[1].intValue() - leafNodes.get(i).lb.data[1].intValue();
                    if (tempY <= 0){ tempY = 1; }
                    int tempX = leafNodes.get(i).ub.data[0].intValue() - leafNodes.get(i).lb.data[0].intValue();
                    if (tempX <= 0){ tempX = 1; }
                    int dummyY = random.nextInt(tempY) + leafNodes.get(i).lb.data[1].intValue();
                    //int dummyX = random.nextInt(tempX)/4 + leafNodes.get(i).lb.data[0].intValue();
                    int dummyX = leafNodes.get(i).lb.data[0].intValue();
//                    int dummyY = random.nextInt(leafNodes.get(i).ub.data[1].intValue() - leafNodes.get(i).lb.data[1].intValue()) + leafNodes.get(i).lb.data[1].intValue();
//                    int dummyX = random.nextInt(leafNodes.get(i).ub.data[0].intValue() - leafNodes.get(i).lb.data[0].intValue())/4 + leafNodes.get(i).lb.data[0].intValue();
                    for (int k = 0; k < maxNum; k++){
                        int predummyX = dummyX;
                        //dummyX = (dummyX + random.nextInt(tempX)/4) % (int)(Math.pow(2,dataLength));
                        dummyX = dummyX + 1;
//                        dummyX = (dummyX + random.nextInt(leafNodes.get(i).ub.data[0].intValue() - leafNodes.get(i).lb.data[0].intValue())/4) % (int)(Math.pow(2,dataLength));
//                        while (dummyX < predummyX){
//                            dummyX = (dummyX + random.nextInt(tempX)/4) % (int)(Math.pow(2,dataLength));
////                            dummyX = (dummyX + random.nextInt(leafNodes.get(i).ub.data[0].intValue() - leafNodes.get(i).lb.data[0].intValue())/4) % (int)(Math.pow(2,dataLength));
//                        }
                        dummyVertices.add(new Vertex(new BigInteger[]{BigInteger.valueOf(dummyX), BigInteger.valueOf(dummyY)}));
                    }

                    VCPointLists.get(i).add(new Point(dummyId, dummyData, dummyVertices));
                }
                idAccumulator += padNum;
            }

            vcs = new VC[vcNum];
            for (int i = 0; i < vcNum; i++) {
                Point low, high;
                low = new Point(leafNodes.get(i).lb.id, leafNodes.get(i).lb.data);
                high = new Point(leafNodes.get(i).ub.id, leafNodes.get(i).ub.data);

                Point[] VCPoints = VCPointLists.get(i).toArray(new Point[VCPointLists.get(i).size()]);
                vcs[i] = new VC(low, high, VCPoints);
            }

            vcsSecrets = shareVCs(vcs, mod);

//            long timeOfIndexCons2 = System.currentTimeMillis() - timeOfIndexConsStart;
//            System.out.println("Index construction time: " + timeOfIndexCons2 + " ms");

        }


        try (Socket socketC1 = new Socket(ipC1, portC1); Socket socketC2 = new Socket(ipC2, portC2);) {
            PrintWriter writerC1 = new PrintWriter(socketC1.getOutputStream());
            BufferedReader readerC1 = new BufferedReader(new InputStreamReader(socketC1.getInputStream()));

            PrintWriter writerC2 = new PrintWriter(socketC2.getOutputStream());
            BufferedReader readerC2 = new BufferedReader(new InputStreamReader(socketC2.getInputStream()));

            Util.writeInt(testType, writerC1);
            Util.writeInt(q, writerC1);
            Util.writeBigInteger(mod, writerC1);
            Util.writeInt(n, writerC1);
            Util.writeInt(d, writerC1);
            if (testType == 0) {
                Util.writePoints(pointsSecretShares[0], writerC1);
            } else if (testType == 1) {
                writerC1.println(parseVCsToJson(vcsSecrets[0]));
            }
            writerC1.println(parseMultiplicationTripleToJson(triples[0]));
            writerC1.println(parseRandomNumberTupleToJson(tuples[0]));
            writerC1.flush();

            Util.writeInt(testType, writerC2);
            Util.writeInt(q, writerC2);
            Util.writeBigInteger(mod, writerC2);
            Util.writeInt(n, writerC2);
            Util.writeInt(d, writerC2);
            if (testType == 0) {
                Util.writePoints(pointsSecretShares[1], writerC2);
            } else if (testType == 1) {
                writerC2.println(parseVCsToJson(vcsSecrets[1]));
            }
            writerC2.println(parseMultiplicationTripleToJson(triples[1]));
            writerC2.println(parseRandomNumberTupleToJson(tuples[1]));
            writerC2.flush();

            Util.writeInt(epoch, writerC1);
            Util.writeInt(epoch, writerC2);

            Point[][] querySecretShares = new Point[2][q];
            for (int k = 0; k < q; k++) {
                querySecretShares[0][k] = new Point(2);
                querySecretShares[1][k] = new Point(2);

                BigInteger[] idSecretShares = randomSplit(queryPoints[k][2], mod);
                querySecretShares[0][k].id = idSecretShares[0];
                querySecretShares[1][k].id = idSecretShares[1];

                for (int j = 0; j < 2; j++) {
                    BigInteger[] pointSecretShares = randomSplit(queryPoints[k][j], mod);

                    querySecretShares[0][k].data[j] = pointSecretShares[0];
                    querySecretShares[1][k].data[j] = pointSecretShares[1];
                }
            }
            for (int i = 0; i < epoch; i++) {
                System.out.print(i + " ");

                Util.writePoints(querySecretShares[0], writerC1);
                Util.writePoints(querySecretShares[1], writerC2);

                int s1 = Util.readInt(readerC1);
                int s2 = Util.readInt(readerC2);

                //System.out.println("s1 = " + s1 + " s2 = " + s2);

                Point[] r1 = Util.readPoints(s1, d, readerC1);
                Point[] r2 = Util.readPoints(s2, d, readerC2);

                if (testType == 0 || testType == 1) {
                    Set<Point> r = new HashSet<>();
                    for (int j = 0; j < s1; j++) {
                        if (r1[j].id.add(r2[j].id).mod(mod).compareTo(BigInteger.ZERO) != 0) {
                            BigInteger[] rdata = new BigInteger[d];
                            for (int k = 0; k < d; k++) {
                                rdata[k] = r1[j].data[k].add(r2[j].data[k]).mod(mod);
                            }
                            r.add(new Point(r1[j].id.add(r2[j].id).mod(mod), rdata));
                        }
                    }
                }

            }
            System.out.println();
            long timeC1 = Util.readLong(readerC1);
            long communicationTimeC1 = Util.readLong(readerC1);
            long computingTimeC1 = Util.readLong(readerC1);
            long timeC2 = Util.readLong(readerC2);
            long communicationTimeC2 = Util.readLong(readerC2);
            long computingTimeC2 = Util.readLong(readerC2);
            System.out.println("Time C1: " + timeC1 + " ms");
            System.out.println("Communication Time C1: " + communicationTimeC1 + " ms");
            System.out.println("Computing Time C1: " + computingTimeC1 + " ms");
            System.out.println("Time C2: " + timeC2 + " ms");
            System.out.println("Communication Time C2: " + communicationTimeC2 + " ms");
            System.out.println("Computing Time C2: " + computingTimeC2 + " ms");
        }
    }

    /**
    * args: role portC1
    * 
    * @param args
    * @throws IOException
    */
    public static void c1(String[] args) throws IOException {

        int index = 1;
        int portC1 = Integer.parseInt(args[index++]);

        ServerSocket serverSocket = new ServerSocket(portC1);

        Socket socketUser = serverSocket.accept();
        PrintWriter writerUser = new PrintWriter(socketUser.getOutputStream());
        BufferedReader readerUser = new BufferedReader(new InputStreamReader(socketUser.getInputStream()));

        Socket socketC2 = serverSocket.accept();
        PrintWriter writerC2 = new PrintWriter(socketC2.getOutputStream());
        BufferedReader readerC2 = new BufferedReader(new InputStreamReader(socketC2.getInputStream()));

        int testType = Util.readInt(readerUser);
        int q = Util.readInt(readerUser);
        BigInteger mod = Util.readBigInteger(readerUser);
        int n = Util.readInt(readerUser);
        int d = Util.readInt(readerUser);
        Point[] points = null;

        VC[] vcs = null;
        if (testType == 0) {
            points = Util.readPoints(n, d, readerUser);
        } else if (testType == 1) {
            vcs = parseJsonToVCs(readerUser.readLine());
        }
        MultiplicationTriple triple = parseJsonToMultiplicationTriple(readerUser.readLine());
        RandomNumberTuple tuple = parseJsonToRandomNumberTuple(readerUser.readLine());

        AdditiveSecretSharing.PartyID partyID = AdditiveSecretSharing.PartyID.C1;

        int epoch = Util.readInt(readerUser);
        long timeSum = 0L;
        long communicationTimeSum = 0l;
        for (int i = 0; i < epoch; i++) {

            Point[] queryPoints = Util.readPoints(q, 2, readerUser);

            Util.writeInt(i, writerC2);
            Util.readInt(readerC2);

            Point[] r1 = null;
            long timePre = System.currentTimeMillis();
            RunningTimeCounter.startRecord(RunningTimeCounter.COMMUNICATION_TIME);
            if (testType == 0) {
                r1 = SkylineQuery.secureBaselineSQ(partyID, points, queryPoints, q, triple, tuple, mod, readerC2, writerC2);
            } else if (testType == 1) {
                r1 = SkylineQuery.secureSQ1(partyID, vcs, queryPoints, q, d, triple, tuple, mod, readerC2, writerC2);
            }

            timeSum += System.currentTimeMillis() - timePre;
            communicationTimeSum += RunningTimeCounter.get(RunningTimeCounter.COMMUNICATION_TIME);

            Util.writeInt(r1.length, writerUser);
            Util.writePoints(r1, writerUser);

        }
        long timeAvg = timeSum / epoch;
        long communicationTimeAvg = communicationTimeSum / epoch;
        long computingTimeAvg = timeAvg - communicationTimeAvg;
        Util.writeLong(timeAvg, writerUser);
        Util.writeLong(communicationTimeAvg, writerUser);
        Util.writeLong(computingTimeAvg, writerUser);

        socketC2.close();
        socketUser.close();
        serverSocket.close();
    }

    /**
    * args: role ipC1 portC1 portC2
    * 
    * @param args
    * @throws IOException
    */
    public static void c2(String[] args) throws IOException {

        int index = 1;
        String ipC1 = args[index++];
        int portC1 = Integer.parseInt(args[index++]);
        int portC2 = Integer.parseInt(args[index++]);

        ServerSocket serverSocket = new ServerSocket(portC2);

        Socket socketUser = serverSocket.accept();
        PrintWriter writerUser = new PrintWriter(socketUser.getOutputStream());
        BufferedReader readerUser = new BufferedReader(new InputStreamReader(socketUser.getInputStream()));

        Socket socketC1 = new Socket(ipC1, portC1);
        PrintWriter writerC1 = new PrintWriter(socketC1.getOutputStream());
        BufferedReader readerC1 = new BufferedReader(new InputStreamReader(socketC1.getInputStream()));

        int testType = Util.readInt(readerUser);
        int q = Util.readInt(readerUser);
        BigInteger mod = Util.readBigInteger(readerUser);
        int n = Util.readInt(readerUser);
        int d = Util.readInt(readerUser);
        Point[] points = null;

        VC[] vcs = null;
        if (testType == 0) {
            points = Util.readPoints(n, d, readerUser);
        } else if (testType == 1) {
            vcs = parseJsonToVCs(readerUser.readLine());
        }
        MultiplicationTriple triple = parseJsonToMultiplicationTriple(readerUser.readLine());
        RandomNumberTuple tuple = parseJsonToRandomNumberTuple(readerUser.readLine());

        AdditiveSecretSharing.PartyID partyID = AdditiveSecretSharing.PartyID.C2;
        int epoch = Util.readInt(readerUser);
        long timeSum = 0L;
        long communicationTimeSum = 0l;
        for (int i = 0; i < epoch; i++) {

            Point[] queryPoints = Util.readPoints(q, 2, readerUser);

            Util.writeInt(i, writerC1);
            Util.readInt(readerC1);

            Point[] r2 = null;
            long timePre = System.currentTimeMillis();
            RunningTimeCounter.startRecord(RunningTimeCounter.COMMUNICATION_TIME);
            // testing function
            if (testType == 0) {
                r2 = SkylineQuery.secureBaselineSQ(partyID, points, queryPoints, q, triple, tuple, mod, readerC1, writerC1);
            } else if (testType == 1) {
                r2 = SkylineQuery.secureSQ1(partyID, vcs, queryPoints, q, d, triple, tuple, mod, readerC1, writerC1);
            }

            timeSum += System.currentTimeMillis() - timePre;
            communicationTimeSum += RunningTimeCounter.get(RunningTimeCounter.COMMUNICATION_TIME);

            Util.writeInt(r2.length, writerUser);
            Util.writePoints(r2, writerUser);

        }
        long timeAvg = timeSum / epoch;
        long communicationTimeAvg = communicationTimeSum / epoch;
        long computingTimeAvg = timeAvg - communicationTimeAvg;
        Util.writeLong(timeAvg, writerUser);
        Util.writeLong(communicationTimeAvg, writerUser);
        Util.writeLong(computingTimeAvg, writerUser);

        socketC1.close();
        socketUser.close();
        serverSocket.close();
    }

    public static Point[] sharePoint(Point point, BigInteger mod) {
        int m = point.data.length;

        BigInteger[] idSecretShares = randomSplit(point.id, mod);
        BigInteger[][] pointDataSecrets = new BigInteger[2][m];
        for (int i = 0; i < m; i++) {
            BigInteger[] secrets = randomSplit(point.data[i], mod);
            pointDataSecrets[0][i] = secrets[0];
            pointDataSecrets[1][i] = secrets[1];
        }

        return new Point[] { new Point(idSecretShares[0], pointDataSecrets[0]),
                new Point(idSecretShares[1], pointDataSecrets[1]) };
    }

    public static Point[] sharePointWithVertex(Point point, BigInteger mod) {
        if (point.vertices == null) {
            return sharePoint(point, mod);
        }
        int m = point.data.length;
        int v = point.vertices.size();

        BigInteger[] idSecretShares = randomSplit(point.id, mod);
        BigInteger[][] pointDataSecrets = new BigInteger[2][m];
        for (int i = 0; i < m; i++) {
            BigInteger[] secrets = randomSplit(point.data[i], mod);
            pointDataSecrets[0][i] = secrets[0];
            pointDataSecrets[1][i] = secrets[1];
        }
        List<Vertex>[] vertexSecrets = new List[2];
        vertexSecrets[0] = new ArrayList<Vertex>();
        vertexSecrets[1] = new ArrayList<Vertex>();
        for (int i = 0; i < v; i++) {
            BigInteger[][] vertexDataSecrets = new BigInteger[2][2];
            for (int j = 0; j < 2; j++){
                BigInteger[] secrets = randomSplit(point.vertices.get(i).data[j], mod);
                vertexDataSecrets[0][j] = secrets[0];
                vertexDataSecrets[1][j] = secrets[1];
            }
            vertexSecrets[0].add(new Vertex(vertexDataSecrets[0]));
            vertexSecrets[1].add(new Vertex(vertexDataSecrets[1]));
        }

        return new Point[] { new Point(idSecretShares[0], pointDataSecrets[0], vertexSecrets[0]),
                new Point(idSecretShares[1], pointDataSecrets[1], vertexSecrets[1]) };
    }

    public static List<Vertex> getVertices(Point point, List<Point> dominatorsList, int id, int dimension, int dataLength, int datanum, int datasetType) {

//        List<Point> pointList = new ArrayList<>(pointListCopy);
        //System.out.println(id);
//        List<Point> dominatorsList = getNonspatialDominators2(pointList, point, id, dimension);
////        System.out.println("dominators size:" + (dominatorsList.size() - 1));
//        if (dominatorsList.isEmpty()) {
//            System.out.println(id);
//            return;
//        }

        double[] dis = new double[dominatorsList.size()];
        for (int i = 0; i < dominatorsList.size(); i++) {
            double dx = point.data[0].doubleValue() - dominatorsList.get(i).data[0].doubleValue();
            double dy = point.data[1].doubleValue() - dominatorsList.get(i).data[1].doubleValue();
            dis[i] = dx * dx + dy * dy;
        }

        int kk = Math.min(dominatorsList.size(), 50);
        int[] indices = findKSmallestIndexes(dis, kk);

//        Integer[] indices = IntStream.range(0, dominatorsList.size()).boxed().sorted(Comparator.comparingDouble(i -> dis[i])).toArray(Integer[]::new);

        int k = Math.min(dominatorsList.size(), 10);
        boolean stop = false;
        List<Vertex> verticesList = new ArrayList<>();
//        long timeOfVoronoiConsStart = System.currentTimeMillis();
        List<Point> subListCopy = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            subListCopy.add(dominatorsList.get(indices[i]));
        }
        int m = 0;
        while (!stop){
            verticesList = VoronoiWithTinfour.getPolygon(point, subListCopy, dataLength);
            if ((k + m) == kk) break;
            Point dp = dominatorsList.get(indices[k + m]);
            BigInteger[] p2Vertex = new BigInteger[verticesList.size()];
            for (int i = 0; i < verticesList.size(); i++) {
                p2Vertex[i] = calDist(point.data, verticesList.get(i).data);
            }
            BigInteger maxDist = Collections.max(Arrays.asList(p2Vertex));
            if (maxDist.multiply(BigInteger.TWO).compareTo(calDist(point.data, dp.data)) < 0) {
                stop = true;
            }else {
                subListCopy.add(dp);
                m += 1;
            }
        }

//        long timeOfVoronoiCons = System.currentTimeMillis() - timeOfVoronoiConsStart;
//        System.out.println("Voronoi construction time: " + timeOfVoronoiCons + " ms");
        //print verticesList
//        for (Vertex vertex : verticesList) {
//            System.out.println(vertex.data[0] + " " + vertex.data[1]);
//        }


//        lock.lock();
//        try {
//            Util.writeTxtFile2(verticesList, id, "src\\main\\resources\\index\\vertices_cc_" + datasetType + "_" + datanum + ".txt");
//        }
//        finally {
//            lock.unlock();
//        }


//        String[][] dominatorsListStr = new String[dominatorsList.size()][2];
//        for (int i = 0; i < dominatorsList.size(); i++) {
//            for (int j = 0; j < 2; j++) {
//                dominatorsListStr[i][j] = dominatorsList.get(i).data[j].toString();
//            }
//        }
//        Util.writeTxtFile(dominatorsListStr, "E:\\secureskylinequery-v2.0-for-large-index\\src\\main\\resources\\index\\points.txt");
//        try{
//            String[] cmd = new String[]{"python", "E:\\secureskylinequery-v2.0-for-large-index\\src\\main\\resources\\getvertex.py"};
//            Process process = Runtime.getRuntime().exec(cmd);
//            process.waitFor();
//        }catch (Exception e){
//            e.printStackTrace();
//        }
//        List<Vertex> verticesList = new ArrayList<>();
//        String[] verticesStr = Util.readTxtFile("E:\\secureskylinequery-v2.0-for-large-index\\src\\main\\resources\\index\\vertices_" + datasetType + "_" + datanum + ".txt");
//        for (int i = 0; i < verticesStr.length; i++) {
//            String[] vertexStr = verticesStr[i].trim().split(" ");
//            BigInteger[] vertexData = new BigInteger[2];
//            for (int j = 0; j < 2; j++) {
//                vertexData[j] = new BigInteger(vertexStr[j]);
//            }
//            verticesList.add(new Vertex(vertexData));
//        }
//
        return verticesList;
    }

    public static AG[][] shareAGs(AG[] ags, BigInteger mod) {
        int num = ags.length;
        int size = ags[0].points.length;
        int m = ags[0].points[0].data.length;

        AG[][] agsSecrets = new AG[2][num];
        for (int i = 0; i < num; i++) {
            agsSecrets[0][i] = new AG(size, m);
            agsSecrets[1][i] = new AG(size, m);

            BigInteger[] labelSecrets = randomSplit(ags[i].label, mod);
            agsSecrets[0][i].label = labelSecrets[0];
            agsSecrets[1][i].label = labelSecrets[1];

            for (int j = 0; j < size; j++) {
                BigInteger[] subLabelSecrets = randomSplit(ags[i].subLabels[j], mod);
                agsSecrets[0][i].subLabels[j] = subLabelSecrets[0];
                agsSecrets[1][i].subLabels[j] = subLabelSecrets[1];

                Point[] pointSecrets = sharePoint(ags[i].points[j], mod);
                agsSecrets[0][i].points[j] = pointSecrets[0];
                agsSecrets[1][i].points[j] = pointSecrets[1];
            }
        }

        return agsSecrets;
    }

    public static String parseAGsToJson(AG[] ags) {
        return JSON.toJSONString(ags);
    }

    public static AG[] parseJsonToAGs(String json) {
        return JSON.parseArray(json, AG.class).toArray(new AG[] {});
    }

    public static VG[][] shareVGs(VG[] vgs, BigInteger mod) {
        int num = vgs.length;
        int size = vgs[0].points.length;
        int m = vgs[0].points[0].data.length;

        VG[][] vgsSecrets = new VG[2][num];
        for (int i = 0; i < num; i++) {
            vgsSecrets[0][i] = new VG(size, m);
            vgsSecrets[1][i] = new VG(size, m);

            Point[] lowSecrets = sharePoint(vgs[i].low, mod);
            vgsSecrets[0][i].low = lowSecrets[0];
            vgsSecrets[1][i].low = lowSecrets[1];

            Point[] highSecrets = sharePoint(vgs[i].high, mod);
            vgsSecrets[0][i].high = highSecrets[0];
            vgsSecrets[1][i].high = highSecrets[1];

            for (int j = 0; j < size; j++) {
                BigInteger[] subLabelSecrets = randomSplit(vgs[i].subLabels[j], mod);
                vgsSecrets[0][i].subLabels[j] = subLabelSecrets[0];
                vgsSecrets[1][i].subLabels[j] = subLabelSecrets[1];

                Point[] pointSecrets = sharePoint(vgs[i].points[j], mod);
                vgsSecrets[0][i].points[j] = pointSecrets[0];
                vgsSecrets[1][i].points[j] = pointSecrets[1];
            }
        }

        return vgsSecrets;
    }

    public static VC[][] shareVCs(VC[] vcs, BigInteger mod) {
        int num = vcs.length;
        int size = vcs[0].points.length;
        int m = vcs[0].points[0].data.length;

        VC[][] vcsSecrets = new VC[2][num];
        for (int i = 0; i < num; i++) {
            vcsSecrets[0][i] = new VC(size, m);
            vcsSecrets[1][i] = new VC(size, m);

            Point[] lowSecrets = sharePoint(vcs[i].low, mod);
            vcsSecrets[0][i].low = lowSecrets[0];
            vcsSecrets[1][i].low = lowSecrets[1];

            Point[] highSecrets = sharePoint(vcs[i].high, mod);
            vcsSecrets[0][i].high = highSecrets[0];
            vcsSecrets[1][i].high = highSecrets[1];

            for (int j = 0; j < size; j++) {

                Point[] pointSecrets = sharePointWithVertex(vcs[i].points[j], mod);
                vcsSecrets[0][i].points[j] = pointSecrets[0];
                vcsSecrets[1][i].points[j] = pointSecrets[1];
            }
        }

        return vcsSecrets;
    }

    public static String parseVGsToJson(VG[] vgs) {
        return JSON.toJSONString(vgs);
    }

    public static VG[] parseJsonToVGs(String json) {
        return JSON.parseArray(json, VG.class).toArray(new VG[] {});
    }

    public static String parseVCsToJson(VC[] vcs) {
        return JSON.toJSONString(vcs);
    }

    public static VC[] parseJsonToVCs(String json) {
        return JSON.parseArray(json, VC.class).toArray(new VC[] {});
    }

    public static void main(String[] args) throws IOException {
        JSON.config(JSONWriter.Feature.LargeObject, true);

        //System.out.println("args: " + Arrays.asList(args));

        int testType = 1; // 0-secskyline  1-ours
        int datasetType = 1; // 0-CORR 1-INDE 2-ANTI 3-HOTE
        int N = 1000; // dataset size
        int dataLength = 8; // the bit length of values
        int dimension = 4; // number of dimensions
        int minVCSize = 30; //the minimum number of entries in a single grid during space partitioning

        String c1 = "c1 8001"; // role portC1
        String c2 = "c2 127.0.0.1 8001 8002"; // role ipC1 portC1 portC2
        String user = "user 127.0.0.1 8001 127.0.0.1 8002 " // role ipC1 portC1 ipC2 portC2
                 + testType + " 10 " + "2024 " + N + " " + dataLength + " " + dimension + " 1 " + datasetType + " " + minVCSize; //testType epoch randomSeed dataNumber dataLength dimension queryNumber dataset minVCSize

        //select a role
        //args = c1.split(" ");
        //args = c2.split(" ");
        args = user.split(" ");

        if (args[0].equals("user"))
            user(args);
        else if (args[0].equals("c1"))
            c1(args);
        else if (args[0].equals("c2"))
            c2(args);
    }
}
