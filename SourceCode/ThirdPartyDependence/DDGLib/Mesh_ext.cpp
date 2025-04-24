//============================================================
// HalfEdge.cpp
// Keenan Crane
// 

#include <map>
#include <cmath>
#include <queue>
#include <float.h>
#include <assert.h>
#include "Mesh.h"
#include "Connection.h"
#include "MeshIO.h"
#include <fstream>
#include <sstream>
#include "qdebug.h"

using namespace std;

namespace DDG
{

   double Mesh::parallelTransport( double phi, HalfEdgeCIter he )
   // given an angle phi relative to the canonical reference frame
   // of he->face, returns the angle parallel transported across he
   // using the Levi-Civita connection, expressed relative to the
   // canonical frame of he->flip->face
   {
      // get (oriented) direction along shared edge
      VertexIter u = he->vertex;
      VertexIter v = he->flip->vertex;
      Vector e = v->position - u->position;
      if( u->index > v->index ) e = -e;

      // compute angle adjustments between canonical frames
      Vector e1, e2; he->face->frame( e1, e2 );
      Vector f1, f2; he->flip->face->frame( f1, f2 );
      double deltaIJ = atan2( dot(e,e2), dot(e,e1) );
      double deltaJI = atan2( dot(e,f2), dot(e,f1) );

      // transport phi
      return ( phi - deltaIJ ) + deltaJI;
   }

   double Mesh::defect( const Cycle& c )
   {
      double theta = 0.;

      for( Cycle::const_iterator he = c.begin(); he != c.end(); he++ )
      {
         theta = parallelTransport( theta, *he );
      }

      while( theta >=  M_PI ) theta -= 2.*M_PI;
      while( theta <  -M_PI ) theta += 2.*M_PI;
      
      return -theta;
   }

   void Mesh::integralCurve( FaceIter initialFace,
                               const Vector& initialPoint,
                               double initialAngle,
                               std::vector<Vector>& curve,
                               std::vector<Vector>& normals,
                               int maxPts )
   {
      static int curveIndex = 0; // unique ID for each curve generated
      FaceIter f = initialFace; // current face
      Vector x = initialPoint; // current point
      double alpha = f->alpha + initialAngle;
      Vector u( cos(alpha), sin(alpha), 0. ); // current direction 

      curve.clear();
      curve.push_back( f->toGlobal( x ));
      normals.push_back( f->normal() );

      for( int i = 0; i < maxPts; i++ )
      {
         // stop if we enter a virtual boundary face
         if( f->isBoundary() ) break;

         // stop if we've already visited this face
         if( f->curveIndex == curveIndex ) break;

         HalfEdgeCIter h[3];
         Vector q[3];

         h[0] = f->he;
         h[1] = f->he->next;
         h[2] = f->he->next->next;
         for( int j = 0; j < 3; j++ )
         {
            q[j] = f->toLocal( h[j]->vertex->position );
         }

         // intersect ray with each edge
         const double eps = 1e-6; // tolerance for intersecting the edge we're coming from
         double tMin = DBL_MAX; // minimum distance to any edge
         HalfEdgeCIter hMin = h[0]; // closest edge
         for( int j = 0; j < 3; j++ )
         {
            int k = (j+1)%3;
            Vector c = q[k]-q[j];
            c = Vector( -c.y, c.x, 0. ).unit();
            double d = -dot(c, q[j]);
            double t = -(dot(c,x)+d)/dot(c,u);

            if( t > eps && t < tMin )
            {
               tMin = t;
               hMin = h[j];
            }
         }

         // stop if intersection test yields dubious results
         if( tMin == DBL_MAX || tMin > 2.*f->circumradius() )
         {
            break;
         }

         // move current point to intersection point
         x += tMin*u;

         curve.push_back( f->toGlobal( x ));
         normals.push_back( f->normal() );

         // get pointer to next triangle
         FaceIter g = hMin->flip->face;

         // rewrite x in local coords of next triangle
         x = g->toLocal( f->toGlobal( x ));

         // rewrite u in local coords of next triangle
         double thetaIJ = hMin->edge->theta;
         if( hMin->vertex->index > hMin->flip->vertex->index ) thetaIJ = -thetaIJ;
         double beta = parallelTransport( atan2(u.y,u.x), hMin ) - thetaIJ;
         u = Vector( cos(beta), sin(beta), 0. );

         // indicate which curve last passed through f
         f->curveIndex = curveIndex;

         // move to next triangle
         f = g;
      }

      curveIndex++; // keep track of which curve we're integrating
   }

   void Mesh::topologyChange( void )
   {
      indexElements();

      buildTreeCotreeDecomposition();
      generatorIndices.resize( nGenerators(), 0. );

      if( connection != NULL )
      {
         delete connection;
         connection = NULL;
      }
      qDebug() << "topologyChange ended!";
   }

   void Mesh::geometryChange( void )
   {
      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         e->updateStar();
      }

      if( connection != NULL )
      {
         delete connection;
         connection = NULL;
      }
   }

   void Mesh::computeFrameAngles( double initialAngle )
   {
      if( transportRoot->constraintAngle >= 0. )
         transportRoot->alpha = transportRoot->constraintAngle;
      else
         transportRoot->alpha = initialAngle;

      for( vector<TransportData>::const_iterator td  = transportOrder.begin();
                                                 td != transportOrder.end();
                                                 td ++ )
      {
         // alphaJ = alphaI + delta - sign*omega
         *(td->alphaJ) = *(td->alphaI) + td->delta - td->sign*(*(td->omega));

         // here we're transporting the angle "alphaI" in triangle i across
         // a shared edge to get the angle "alphaJ" in neighboring triangle j;
         // "delta" compensates for the difference between the local
         // reference frames, "omega" is the angle of the connection, and
         // "sign" accounts for the fact that the direction of transport
         // might not be consistent with the orientation of the shared
         // dual edge.
      }
   }

   bool Mesh::hasBoundary( void ) const
   {
      for( HalfEdgeCIter he = halfedges.begin();
                         he != halfedges.end();
                         he++ )
      {
         if( he->onBoundary )
         {
            return true;
         }
      }
      return false;
   }

   bool inPrimalSpanningTree( HalfEdgeCIter he )
   {
      VertexCIter v = he->vertex;
      VertexCIter w = he->flip->vertex;

      return v->parent == w || w->parent == v;
   }

   bool inDualSpanningTree( HalfEdgeCIter he )
   {
      FaceCIter f = he->face;
      FaceCIter g = he->flip->face;
     
      return f->parent == g || g->parent == f;
   }

   void Mesh::buildPrimalSpanningTree( void )
   {
      VertexIter root = vertices.begin();
      while( root->onBoundary() ) root++;

      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->parent = v;
      }

      queue<VertexIter> Q;
      Q.push( root );
      while( !Q.empty() )
      {
         VertexIter v = Q.front(); Q.pop();

         HalfEdgeIter he = v->he;
         do
         {
            VertexIter w = he->flip->vertex;

            if( w->parent == w &&
                w != root &&
                !w->onBoundary() )
            {
               w->parent = v;
               Q.push( w );
            }

            he = he->flip->next;
         }
         while( he != v->he );
      }
   }

   void Mesh::buildDualSpanningCoTree( void )
   {
      FaceIter root = faces.begin();
      while( root->isBoundary() ) root++;

      for( FaceIter f = faces.begin(); f != faces.end(); f++ )
      {
         f->parent = f;
      }

      queue<FaceIter> Q;
      Q.push( root );
      while( !Q.empty() )
      {
         FaceIter f = Q.front(); Q.pop();

         HalfEdgeIter he = f->he;
         do
         {
            FaceIter g = he->flip->face;

            if( g->parent == g &&
                g != root &&
                !inPrimalSpanningTree( he ) &&
                !g->isBoundary() )
            {
               g->parent = f;
               Q.push( g );
            }

            he = he->next;
         }
         while( he != f->he );
      }
   }

   void Mesh::buildTreeCotreeDecomposition( void )
   {
      buildPrimalSpanningTree();
      buildDualSpanningCoTree();
   }

   HalfEdgeIter sharedHalfEdge( VertexIter& v, VertexIter& w )
   {
      HalfEdgeIter he = v->he;
      do
      {
         if( he->flip->vertex == w )
         {
            return he;
         }

         he = he->flip->next;
      }
      while( he != v->he );
      
      assert( 0 );
   }

   HalfEdgeIter sharedHalfEdge( FaceIter& f, FaceIter& g )
   {
      HalfEdgeIter he = f->he;
      do
      {
         if( he->flip->face == g )
         {
            return he;
         }

         he = he->next;
      }
      while( he != f->he );
      
      assert( 0 );
   }

   int Mesh::nGenerators( void ) const
   {
      int n = 0;

      for( EdgeCIter e = edges.begin(); e != edges.end(); e++ )
      {
         if( e->onBoundary() ) continue;

         if( !inPrimalSpanningTree( e->he ) &&
             !inDualSpanningTree( e->he ))
         {
            n++;
         }
      }

      return n;
   }

   bool Mesh::isDualBoundaryLoop( const Cycle& cycle )
   {
      if( cycle.size() == 0 ) return false;

      return cycle[0]->vertex->onBoundary() ||
             cycle[0]->flip->vertex->onBoundary();
   }

   void Mesh::appendDualGenerators( vector<Cycle>& cycles )
   {
      for( EdgeIter e = edges.begin(); e != edges.end(); e++ )
      {
         if( e->onBoundary() ) continue;

         if( !inPrimalSpanningTree( e->he ) &&
             !inDualSpanningTree( e->he ))
         {
            Cycle g, c1, c2;
            FaceIter f;

            g.push_back( e->he );

            f = e->he->flip->face;
            while( f != f->parent )
            {
               c1.push_back( sharedHalfEdge( f, f->parent ));
               f = f->parent;
            }

            f = e->he->face;
            while( f != f->parent )
            {
               c2.push_back( sharedHalfEdge( f, f->parent ));
               f = f->parent;
            }

            int m = c1.size()-1;
            int n = c2.size()-1;
            while( c1[m] == c2[n] ) { m--; n--; }
            for( int i = 0; i <= m; i++ ) g.push_back( c1[i] );
            for( int i = n; i >= 0; i-- ) g.push_back( c2[i]->flip );
            
            // make sure that boundary loops wind around the boundary in a consistent direction
            if( isDualBoundaryLoop( g ))
            {
               if( g[0]->next->vertex->onBoundary() )
               {
                  unsigned int n = g.size();
                  for( unsigned int i = 0; i < n; i++ )
                  {
                     g[i] = g[i]->flip;
                  }

                  for( unsigned int i = 0; i < n/2; i++ )
                  {
                     swap( g[i], g[n-1-i] );
                  }
               }
            }

            cycles.push_back( g );
         }
      }
   }

   void Mesh::appendDirectionalConstraints( vector<Cycle>& cycles, vector<double>& holonomies )
   {
      // first point all faces to themselves to indicate that they have not yet
      // been added to the tree; meanwhile look for a constrained face to serve
      // as the root for our constraint tree (if there aren't any constrained
      // faces, an arbitrary face will work just fine)
      transportRoot = faces.begin();
      for( FaceIter f = faces.begin(); f != faces.end(); f++ )
      {
         f->cParent = f;
         if( f->constraintAngle >= 0. )
         {
            transportRoot = f;
         }
      }

      // do a breadth first search on the dual edges and cache the traversal
      // order and angles (we'll need this tree later to construct a global
      // frame starting with a known direction)
      queue<FaceIter> Q;
      Q.push( transportRoot );
      while( !Q.empty() )
      {
         FaceIter f = Q.front(); Q.pop();

         // visit neighboring faces
         HalfEdgeIter he = f->he;
         do
         {
            FaceIter g = he->flip->face;
            if( g->cParent == g &&
                g != transportRoot &&
                !g->isBoundary() )
            {
               // point the current neighbor to its parent in the
               // traversal and enqueue it
               g->cParent = f;
               Q.push( g );

               // also, cache transport via Levi-Civita across this
               // edge (can then transport via the trivial connection
               // by simply subtracting connection angles once we
               // compute them)
               TransportData td;
               td.delta = parallelTransport( 0., he );
               td.sign = he->vertex->index > he->flip->vertex->index ? -1. : 1.;
               td.omega = &(he->edge->theta);
               td.alphaI = &(he->face->alpha);
               td.alphaJ = &(he->flip->face->alpha);
               transportOrder.push_back( td );

               // check if this face is constrained; if so, we need to add a constraint
               if( g->constraintAngle >= 0. )
               {
                  Cycle constraint;
                  double alpha = g->constraintAngle;

                  // follow this face back up the tree to the most
                  // recent constrained ancestor, adding all the halfedges
                  // in between to a new constraint; also compute the difference
                  // between the two constrained directions relative to transport
                  // via the Levi-Civita connection
                  FaceIter h = g;
                  do
                  {
                     HalfEdgeIter he = sharedHalfEdge( h, h->cParent );
                     alpha = parallelTransport( alpha, he );
                     constraint.push_back( he );
                     h = h->cParent;
                  } while( h->constraintAngle < 0. );

                  // add this new constraint to the set of all constraints
                  cycles.push_back( constraint );
                  double gamma = h->constraintAngle;
                  Vector u1( cos(gamma), sin(gamma), 0. );
                  Vector u2( cos(alpha), sin(alpha), 0. );
                  holonomies.push_back( acos(dot(u1,u2)) );
               }
            }
            he = he->next;
         } while( he != f->he );
      }
   }

   double tipAngle( const Vector& x, const Vector& a, const Vector& b )
   // returns the angle between (a-x) and (b-x)
   {
      Vector u = ( a - x ).unit();
      Vector v = ( b - x ).unit();

      return atan2( cross(u, v).norm(), dot(u,v) );
   }

   double Mesh::boundaryLoopCurvature( const Cycle& cycle )
   {
      double totalK = 0.;

      // get a halfedge of the "virtual" face bounded by the current cycle
      VertexCIter v0 = cycle[0]->flip->next->vertex;
      HalfEdgeIter he0 = v0->he;
      do
      {
         he0 = he0->flip->next;
      }
      while( !he0->onBoundary );

      // compute a "virtual" vertex in the middle of this loop
      Vector c( 0., 0., 0. );
      HalfEdgeCIter he = he0;
      int boundaryLength = 0;
      do
      {
         c += he->vertex->position;
         boundaryLength++;
         he = he->next;
      }
      while( he != he0 );
      c /= (double) boundaryLength;

      // compute the curvature around the center vertex
      double K = 2.*M_PI;
      he = he0;
      do
      {
         Vector a = he->vertex->position;
         Vector b = he->next->vertex->position;
         K -= tipAngle( c, a, b );
         he = he->next;
      }
      while( he != he0 );
      totalK += K;

      // add the curvature around each of the boundary vertices, using
      // the following labels:
      //    c - virtual center vertex of boundary loop (computed above)
      //    d - current boundary vertex (we walk around the 1-ring of this vertex)
      //    a,b - consecutive interior vertices in 1-ring of d
      //    e,f - boundary vertices adjacent to d
      he = he0;
      do
      {
         VertexCIter v = he->vertex;
         Vector d = v->position;

         K = 2.*M_PI;

         HalfEdgeCIter he2 = v->he;
         do
         {
            if( he2->onBoundary )
            {
               Vector f = he2->next->vertex->position;
               K -= tipAngle( d, f, c );
            }
            else
            {
               Vector a = he2->next->vertex->position;
               Vector b = he2->next->next->vertex->position;
               K -= tipAngle( d, a, b );

               if( he2->flip->onBoundary )
               {
                  Vector e = he2->flip->vertex->position;
                  K -= tipAngle( d, c, e );
               }
            }

            he2 = he2->flip->next;
         }
         while( he2 != v->he );

         totalK += K;

         he = he->next;
      }
      while( he != he0 );

      return totalK;
   }

   void Mesh::computeTrivialConnection( void )
   {
      if( connection == NULL )
      {
         connection = new Connection( *this );
      }

      connection->update();
      computeFrameAngles( fieldAngle );
   }

   string getFilenameExtension( const string& filename )
   {
      int lastDot = filename.find_last_of( '.' );
      int extLength = filename.length() - lastDot - 1;

      string extension = filename.substr( lastDot+1, extLength );

      transform( extension.begin(), extension.end(), extension.begin(), ::tolower );

      return extension;
   }

   //int Mesh::eulerCharacteristic( void ) const
   //{
   //   return vertices.size() - edges.size() + faces.size();
   //}


   void Mesh::writeEOBJ(std::string path)
   {
       ofstream out(path);
       out.precision(10);

       int currentIndex = 1;
       map< VertexCIter, int > vertexIndex;
       map< FaceCIter, int > faceIndex;
       const vector<Vertex>& vertices(vertices);
       const vector<Face>& faces(faces);

       for (VertexCIter i = vertices.begin(); i != vertices.end(); i++)
       {
           out << "v " << i->position[0] << " "
               << i->position[1] << " "
               << i->position[2] << endl;

           vertexIndex[i] = currentIndex;
           currentIndex++;
       }

       currentIndex = 1;
       for (FaceCIter i = faces.begin(); i != faces.end(); i++)
       {
           HalfEdgeIter he = i->he;

           // don't write boundary faces
           if (he->onBoundary)
           {
               continue;
           }

           faceIndex[i] = currentIndex;
           currentIndex++;

           out << "f ";

           do
           {
               out << vertexIndex[he->vertex] << " ";
               he = he->next;
           } while (he != i->he);

           out << endl;
       }

       for (FaceCIter i = faces.begin(); i != faces.end(); i++)
       {
           HalfEdgeIter he = i->he;

           // don't write vectors on boundary faces
           if (he->onBoundary)
           {
               continue;
           }

           out << "# attrs f " << faceIndex[i] << " ";

           double alpha = i->alpha;
           Vector w(cos(alpha), sin(alpha), 0.);
           Vector u = i->toGlobal(w);

           out << u.x << " " << u.y << " " << u.z;

           out << endl;
       }
   }
}


