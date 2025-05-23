#include <vector>
using namespace std;

#include "Vertex.h"
#include "Mesh.h"
#include "HalfEdge.h"
#include "Utility.h"

namespace DDG
{
   Vector Vertex::normal( void ) const
   // returns the vertex normal
   {
      Vector N;

      HalfEdgeCIter h = he;
      do
      {
         if( !h->onBoundary )
         {
            N += h->face->normal();
         }

         h = h->flip->next;
      }
      while( h != he );

      return N.unit();
   }

   double Vertex::dualArea( void ) const
   {
      double A = 0.;

      HalfEdgeCIter h = he;
      do
      {
         if( !h->onBoundary )
         {
            A += h->face->area();
         }

         h = h->flip->next;
      }
      while( h != he );

      return A/3.;
   }

   vector<HalfEdge> isolated; // all isolated vertices point to isolated.begin()

   bool Vertex::onBoundary( void ) const
   {
       HalfEdgeCIter h = he;
       if (he->onBoundary) {
           return true;
       }
       do
       {
           h = h->next->next->flip;
       } while (h != he && (!h->onBoundary));
       if (h == he) {
           return false;
       }
       else {
           return true;
       }
      //return he->onBoundary;
   }

   bool Vertex::isIsolated( void ) const
   // returns true if the vertex is not contained in any face or edge; false otherwise
   {
      return isolated.size()!=0 && he == isolated.begin();
   }

   int Vertex :: valence( void ) const
   // returns the number of incident faces
   {
      int n = 0;

      HalfEdgeCIter h = he;
      do
      {
         n++;
         h = h->flip->next;
      }
      while( h != he );

      return n;
   }

   double Vertex :: angleSum( void ) const
   {
      double sum = 0.;

      HalfEdgeCIter h = he;
      do
      {
         sum += h->next->angle();
         h = h->flip->next;
      }
      while( h != he );

      return sum;
   }

   double Vertex::defect(void) const
   {
       double sum = 0.;

       // iterate over incident triangles
       HalfEdgeIter he = this->he;
       do
       {
           // grab vertices
           Vector p1 = he->vertex->position;
           Vector p2 = he->next->vertex->position;
           Vector p3 = he->next->next->vertex->position;

           // subtract incident angle from sum
           Vector u1 = (p2 - p1);
           Vector u2 = (p3 - p1);
           sum += atan2(cross(u1, u2).norm(), dot(u1, u2));

           he = he->flip->next;
   }
       while (he != this->he);

   return 2. * M_PI - sum;
   }

   Vector Vertex :: fieldVector( double k, double n ) const
   {
#ifdef SP_FLAT_TORUS
      Complex p = he->origTexCoord;
      Complex X = he->next->origTexCoord - p;
      X /= X.norm();
      double r = directionField.norm();
      double theta = directionField.arg() / k;
      double phi = theta + (2.*M_PI/k)*n;
      Complex Y = Complex( cos(phi), sin(phi) ) * X;
      return r * Vector( Y.re, Y.im, 0. );
#else
      Vector p = position;
      Vector N = normal();
      Vector X = he->flip->vertex->position - p;
      X -= dot(X,N)*N;
      X.normalize();
      Vector JX = cross( N, X );
      double r = directionField.norm();
      double theta = directionField.arg() / k;
      double phi = theta + (2.*M_PI/k)*n + M_PI/2.;
      return r * ( cos(phi)*X + sin(phi)*JX );
#endif
   }

   Complex Vertex :: canonicalVector( void ) const
   {
      double r = directionField.norm();
      double theta = directionField.arg() / 2.;
      return r * Complex( cos(theta), sin(theta) );
   }

   Complex Vertex :: principalDirection( void ) const
   // This is literally the only routine in the entire algorithm
   // that does anything involving principal curvatures/curvature estimates.
   // No thresholds, no parameters...
   {
      Complex X( 0., 0. );

      HalfEdgeCIter h = he;
      do
      {
         double l = h->edge->length();
         double alpha = h->edge->dihedralAngle();
         double theta = h->angularCoordinate;
         Complex r( cos(2.*theta), sin(2.*theta) );

         X += l * alpha * r;

         h = h->flip->next;
      }
      while( h != he );

      return -X / 4.;
   }
}

