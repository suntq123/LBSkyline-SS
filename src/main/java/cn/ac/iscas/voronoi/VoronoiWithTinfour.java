package cn.ac.iscas.voronoi;

//import java.awt.geom.Rectangle2D;
import cn.ac.iscas.skylinequery.SkylineQuery;
import org.tinfour.standard.IncrementalTin;
import org.tinfour.common.INeighborhoodPointsCollector;
import org.tinfour.voronoi.BoundedVoronoiBuildOptions;
import org.tinfour.voronoi.BoundedVoronoiDiagram;
import org.tinfour.common.Vertex;
import org.tinfour.voronoi.ThiessenPolygon;
import org.tinfour.common.IIncrementalTin;
import org.tinfour.common.IQuadEdge;

import java.awt.geom.Rectangle2D;
import java.math.BigInteger;
import java.util.*;

public class VoronoiWithTinfour {

    public static List<SkylineQuery.Vertex> getPolygon(SkylineQuery.Point point, List<SkylineQuery.Point> dominatorsList, int dataLength) {
        List<SkylineQuery.Vertex> verticesList = new ArrayList<>();

        dominatorsList.add(0, point);
        List<Vertex> points = new ArrayList<>();
        for (int i = 0; i < dominatorsList.size(); i++) {
            Vertex vertex = new Vertex(dominatorsList.get(i).data[0].doubleValue(), dominatorsList.get(i).data[1].doubleValue(), 0);
            points.add(vertex);
        }

        // int k = 500;
        // if(points.size() > k){
        //     Vertex p = points.get(0);
        //     points.remove(0);
        //     List<Vertex> newPoints = findNearestPoints(p, points, k);
        //     points.clear();
        //     points.addAll(newPoints);
        //     points.add(0, p);
        // }

//        int k =6;
//        Vertex p = points.get(0);
//        IIncrementalTin tin = new IncrementalTin();
//        tin.add(points, null);
//
//        //List<Vertex> neighbors = tin.getAdjacentVertices(p);
//        INeighborhoodPointsCollector npc = tin.getNeighborhoodPointsCollector();
//        List<Vertex> neighbors = npc.collectNeighboringVertices(p.getX(), p.getY(), 2, k);
//        neighbors.add(0, p);
//        points.clear();
//        points.addAll(neighbors);

        double boundUp = Math.pow(2, dataLength);
        //left
        List<Vertex> pointsl = new ArrayList<>();
        for (int i = 0; i < points.size(); i++) {
            pointsl.add(new Vertex(0 - (points.get(i).getX() - 0), points.get(i).getY(), 0));
        }
        //right
        List<Vertex> pointsr = new ArrayList<>();
        for (int i = 0; i < points.size(); i++) {
            pointsr.add(new Vertex(boundUp + (boundUp - points.get(i).getX()), points.get(i).getY(), 0));
        }
        //down
        List<Vertex> pointsd = new ArrayList<>();
        for (int i = 0; i < points.size(); i++) {
            pointsd.add(new Vertex(points.get(i).getX(), 0 - (points.get(i).getY() - 0), 0));
        }
        //up
        List<Vertex> pointsu = new ArrayList<>();
        for (int i = 0; i < points.size(); i++) {
            pointsu.add(new Vertex(points.get(i).getX(), boundUp + (boundUp - points.get(i).getY()), 0));
        }

        points.addAll(pointsl);
        points.addAll(pointsr);
        points.addAll(pointsd);
        points.addAll(pointsu);

//        int k = 10;
//        Vertex p = points.get(0);
//        IIncrementalTin tin = new IncrementalTin();
//        tin.add(points, null);
//
//        //List<Vertex> neighbors = tin.getAdjacentVertices(p);
//        INeighborhoodPointsCollector npc = tin.getNeighborhoodPointsCollector();
//        List<Vertex> neighbors = npc.collectNeighboringVertices(p.getX(), p.getY(), 2, k);
//        neighbors.add(0, p);
//        points.clear();
//        points.addAll(neighbors);

        BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
        //options.setBounds(new Rectangle2D.Double(0, boundUp, boundUp + 0.5, boundUp + 0.5));

        BoundedVoronoiDiagram voronoi = new BoundedVoronoiDiagram(points, options);

        ThiessenPolygon polygon = voronoi.getContainingPolygon(points.get(0).getX(), points.get(0).getY());
        List<IQuadEdge> edges = polygon.getEdges();
        for (int i = edges.size()-1; i >= 0; i--){
            IQuadEdge e = edges.get(i);
            Vertex v1 = e.getA();
            v1 = fixOutofBoundaryPoint(v1, dataLength);
            verticesList.add(new SkylineQuery.Vertex(new BigInteger[]{BigInteger.valueOf(Math.round(v1.getX())), BigInteger.valueOf(Math.round(v1.getY()))}));
        }

        return  verticesList;
    }

    public static Vertex fixOutofBoundaryPoint(Vertex point, int dataLength) {
        double x, y;
        x = point.getX();
        y = point.getY();
        double boundUp = Math.pow(2, dataLength);
        if (x < 0) {
            x = 0;
        } else if (x > boundUp) {
            x = boundUp;
        }
        if (y < 0) {
            y = 0;
        } else if (y > boundUp) {
            y = boundUp;
        }
        return new Vertex(x, y, 0);
    }

    private static class PointWithDistance {
        Vertex point;
        double distanceSq;

        PointWithDistance(Vertex point, Vertex p) {
            this.point = point;
            double dx = point.getX() - p.getX();
            double dy = point.getY() - p.getY();
            this.distanceSq = dx * dx + dy * dy; // 计算距离平方
        }
    }

    public static List<Vertex> findNearestPoints(Vertex p, List<Vertex> points, int k) {
        if (k <= 0 || points == null || points.isEmpty()) {
            return Collections.emptyList();
        }

        PriorityQueue<PointWithDistance> maxHeap = new PriorityQueue<>(
                Comparator.comparingDouble((PointWithDistance pwd) -> pwd.distanceSq).reversed()
        );

        for (Vertex point : points) {
            PointWithDistance pwd = new PointWithDistance(point, p);
            if (maxHeap.size() < k) {
                maxHeap.offer(pwd);
            } else {
                if (pwd.distanceSq < maxHeap.peek().distanceSq) {
                    maxHeap.poll();
                    maxHeap.offer(pwd);
                }
            }
        }

        List<PointWithDistance> sortedByDistance = new ArrayList<>(maxHeap);
        sortedByDistance.sort(Comparator.comparingDouble(pwd -> pwd.distanceSq));

        List<Vertex> result = new ArrayList<>(k);
        for (PointWithDistance pwd : sortedByDistance) {
            result.add(pwd.point);
        }

        return result;
    }


    public static void main(String[] args) {
        List<Vertex> points = new ArrayList<>();
        points.add(new Vertex(0, 0, 0));
        points.add(new Vertex(0, 1, 0));
        points.add(new Vertex(0, 2, 0));
        points.add(new Vertex(1, 0, 0));
        points.add(new Vertex(1, 1, 0));
        points.add(new Vertex(1, 2, 0));
        points.add(new Vertex(2, 0, 0));
        points.add(new Vertex(2, 1, 0));
        points.add(new Vertex(2, 2, 0));

//        IncrementalTin tin = new IncrementalTin();
//        tin.add(points, null);

        BoundedVoronoiBuildOptions options = new BoundedVoronoiBuildOptions();
//        options.setBounds();
//        options.setBounds(0, 0, 7, 7); // (xMin, yMin, xMax, yMax)

        BoundedVoronoiDiagram voronoi = new BoundedVoronoiDiagram(points, options);

        ThiessenPolygon polygon = voronoi.getContainingPolygon(points.get(4).getX(), points.get(4).getY());
        System.out.println("Site: (" + points.get(4).getX() + ", " + points.get(4).getY() + ")");
        System.out.println("Vertices:");
        List<IQuadEdge> edges = polygon.getEdges();
        for (int i = edges.size()-2; i >= 0; i--){
            IQuadEdge e = edges.get(i);
            Vertex v1 = e.getA();
            System.out.println(" (" + v1.getX() + ", " + v1.getY() + ")");
        }
//        for (IQuadEdge e : edges) {
//            Vertex v1 = e.getA();
//            //Vertex v2 = e.getB();
//            System.out.println(" (" + v1.getX() + ", " + v1.getY() + ")");
//        }

//        List<Vertex> points = new ArrayList<>();
//        Random rand = new Random();
//        for (int i = 0; i < 100000; i++) {
//            points.add(new Vertex(rand.nextDouble() * 1000, rand.nextDouble() * 1000, 0));
//        }
//
//        Vertex p = new Vertex(50.0, 50.0, 0);
//
//        long startTime = System.currentTimeMillis();
//        List<Vertex> nearestPoints = findNearestPoints(p, points, 10);
//        long endTime = System.currentTimeMillis();
//
//        System.out.println("time: " + (endTime - startTime) + "ms");
//        System.out.println("k10:");
//        for (int i = 0; i < 10; i++) {
//            System.out.println(nearestPoints.get(i));
//        }
    }
}
