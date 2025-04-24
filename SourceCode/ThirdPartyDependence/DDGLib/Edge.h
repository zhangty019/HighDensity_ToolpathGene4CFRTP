// -----------------------------------------------------------------------------
// libDDG -- Edge.h
// -----------------------------------------------------------------------------
//
// Edge stores attributes associated with a mesh edge.  The iterator he points
// to one of its two associated halfedges.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_EDGE_H
#define DDG_EDGE_H

#include "Types.h"
#include "HalfEdge.h"
#include "Vector.h"
#include "Vertex.h"

namespace DDG
{
   class Edge
   {
      public:
         HalfEdgeIter he;
         // points to one of the two halfedges associated with this edge

         double length( void ) const;
         // returns the edge length

         double dihedralAngle( void ) const;
         // returns signed dihedral angle

         int index;
         // unique ID in the range [0,nE-1]

         double omega;
         // 1-form guiding parameterization

         bool crossesSheets;
         // whether the target coordinate is conjugated

         // add for paper
         double theta; // connection between incident faces
         double star;  // Hodge star on primal 1-forms
         bool onBoundary(void) const
         {
             return he->onBoundary || he->flip->onBoundary;
         }
         void Edge::updateStar(void)
         {
             double sum = 0.;

             // compute the dual/primal length ratio by adding up the
             // cotangents of the two angles opposing the current edge
             HalfEdgeCIter h = he;
             do
             {
                 Vector& a(h->vertex->position);
                 Vector& b(h->next->vertex->position);
                 Vector& c(h->next->next->vertex->position);

                 // compute the cotangent of the angle of the current
                 // triangle that opposes the current edge
                 Vector u = a - c;
                 Vector v = b - c;
                 double cotTheta = dot(u, v) / cross(u, v).norm();

                 sum += .5 * cotTheta;

                 h = h->flip;
             } while (h != he);

             star = std::max(sum, 0.);
         }
   };
}

#endif
