///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.cpp
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>
#include <unistd.h>
#include <sys/stat.h>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024
using namespace std;
const static int  NOT_FOUND = -1;

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
    : centreOfGravity(0.0,0.0,0.0)
    { // DirectedEdgeSurface()
    // force arrays to size 0
    vertices.resize(0);
    normals.resize(0);
	firstDirectedEdge.resize(0);
	faceVertices.resize(0);
	otherHalf.resize(0);
    } // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (vertexID != vertices.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
			
			// and add it to the vertices
			vertices.push_back(newVertex);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (normalID != normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
			normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (FDEID != firstDirectedEdge.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
			
			// and add it to the vertices
			firstDirectedEdge.push_back(newFDE);
			} // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (faceID != faceVertices.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
			faceVertices.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
			if (otherHalfID != otherHalf.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new face vertex (3 times)
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;
			otherHalf.push_back(newOtherHalf);
			} // other half
        } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
        { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];
        
        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();         
            
            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code
    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
	geometryStream << "# Created for Leeds COMP 5821M Autumn 2020" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "#" << std::endl; 
	geometryStream << "# Surface vertices=" << vertices.size() << " faces=" << faceVertices.size()/3 << std::endl; 
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << vertices[vertex] << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < firstDirectedEdge.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << firstDirectedEdge[vertex] << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < faceVertices.size(); )
        { // per face
        geometryStream << "Face " << face/3 << " ";
        
        // read in three vertices
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++] << " ";
        geometryStream << faceVertices[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
	for (unsigned int dirEdge = 0; dirEdge < otherHalf.size(); dirEdge++)
		geometryStream << "OtherHalf " << dirEdge << " " << otherHalf[dirEdge] << std::endl;
    } // WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= objectSize;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-centreOfGravity.x, -centreOfGravity.y, -centreOfGravity.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
	for (unsigned int face = 0; face < faceVertices.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
			Cartesian3 pq = vertices[faceVertices[face+1]] - vertices[faceVertices[face]];
			Cartesian3 pr = vertices[faceVertices[face+2]] - vertices[faceVertices[face]];

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
					normals[faceVertices[vertex]].x * scale,
					normals[faceVertices[vertex]].y * scale,
					normals[faceVertices[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
				vertices[faceVertices[vertex]].x,
				vertices[faceVertices[vertex]].y,
				vertices[faceVertices[vertex]].z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
	for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glTranslatef(vertices[vertex].x, vertices[vertex].y, vertices[vertex].z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
    } // Render()

// the outline of my strategy: ( I use method 3 taught in our last lab )
// 1.generate a new vertex from every old edge;
//   append the new vertex coordinates to the vector vertices;
//   append the new faces vertex indexes to the vector faceVertices, here, we just append
//   the new faces which have 2 new vertexes and 1 old vertex, in other words, they are not
//   the center new triangles;
//   store center new faces that are formed by three new vertexes in a new vector called
//   newCenterFaces;
// 2.update all of the old vertexes coordinates;
// 3.update all old face vertex indexes;
// 4.generate new first directed edge data;
// 5.generate new normals, which is nor required but I tried.

//          c
//       /    \
//     /        \
//    a - - m - - b
//     \        /
//      \      /
//        \  /
//         d
/***
 * the core function of processing loop subdivision
 */
void DirectedEdgeSurface::loopSubDivision()
{
    //as we need to know which vertexes are old ones, so we
    //store the number of old vertexes. like a milestone.
    int oldVerticesNum = vertices.size();
    //so does the number of old faces.
    int oldFacesNum = faceVertices.size() / 3;
    //so does the number of old face vertex.
    int oldFaceVerticesLength = faceVertices.size();
    //we store all faces fromed by 3 new vertexes, the center face, in the newCenterFaces.
    std::vector<int > newCenterFaces;
    //when we traverse the half edges to create new vertexes,
    //we may find a pair of half edges, but they just have the same new vertex,
    //so we need to know if an edge has created a new vertex, its half edge does not
    //need to create a new vertex. Hence, we use edgesHaveNewVertex to tag if an edge already
    //has a new vertex. Besides, we also store the index to an edge and its half.
    std::vector<int> edgesHaveNewVertex;
    edgesHaveNewVertex.resize(otherHalf.size(),NOT_FOUND);

    //via this loop we can calculate all new vertexes.
    for (int face_id = 0; face_id < oldFacesNum; face_id++) {
        for (int i = 0; i < 3; ++i) {
            //inside a face, we need to check and calculate the new vertexes of 3 edges.
            //just like a circle.
            //here we generate the new vertex on the edge a-b, so
            //we also find vertex d which is the opposite one of c vertex.
            //a and b have 3/8 weight; c and d have 1/8 weight.
            int a_index = face_id * 3 + i % 3;
            int b_index = face_id * 3 + (i + 1) % 3;
            int c_index = face_id * 3 + (i + 2) % 3;

            int b_other_half_edge = otherHalf[b_index];
            //if the edge a-b has already created a new vertex, wo need to jump over it.
            if (edgesHaveNewVertex[b_index] == NOT_FOUND && edgesHaveNewVertex[b_other_half_edge] == NOT_FOUND)
            {
                //d vertex index is always the next one of b's other half edge.
                int d_index = getNextEdge(b_other_half_edge);
                Cartesian3 new_vertex_on_d_coor;
                //calculate the new vertex's coordinates.
                new_vertex_on_d_coor = 0.375f * (getVertexCoor(a_index) + getVertexCoor(b_index))
                                       + 0.125 * (getVertexCoor(c_index) + getVertexCoor(d_index));
                //push it directly to the vector vertices. as we know that which are old ones.
                vertices.push_back(new_vertex_on_d_coor);
                int newVertexIndex = vertices.size() - 1;
                //remember which edges have already created the new vertex.
                edgesHaveNewVertex[b_index] = newVertexIndex;
                edgesHaveNewVertex[b_other_half_edge] = newVertexIndex;
                //every new vertex is need to store in newCenterFaces, and later we will use it to update old faces.
                newCenterFaces.push_back(newVertexIndex);
                //2 new vertexes and 1 old vertex we need to append them to faceVertices.
                faceVertices.push_back(newVertexIndex);

            }//if (edgesHaveNewVertex[b_index] != NOT_FOUND && edgesHaveNewVertex[b_other_half_edge] != NOT_FOUND)
            else
            {
                //the two if statements below are just a guarantee
                int newVertexIndex = edgesHaveNewVertex[b_index];
                if (newVertexIndex == NOT_FOUND)
                {
                    newVertexIndex = edgesHaveNewVertex[b_other_half_edge];
                    edgesHaveNewVertex[b_index] = newVertexIndex;
                }

                if (edgesHaveNewVertex[b_other_half_edge] == NOT_FOUND)
                {
                    edgesHaveNewVertex[b_other_half_edge] = newVertexIndex;
                }
                //every new vertex is need to store in newCenterFaces, and later we will use it to update old faces.
                newCenterFaces.push_back(newVertexIndex);
                //2 new vertexes and 1 old vertex we need to append them to faceVertices.
                faceVertices.push_back(newVertexIndex);
            }
            if (faceVertices.size() % 3 == 0)
            {
                //as every two triangles have one connected vertex,
                //so when we finish form a triangle we need to push the last vertex again.
                faceVertices.push_back(faceVertices[faceVertices.size()-1]);
            }
            //2 new vertexes and 1 old vertex we need to append them to faceVertices.
            faceVertices.push_back(faceVertices[b_index]);
            if (i==2)
            {
                //without this, the last new triangle will miss the last vertex.
                faceVertices.push_back(edgesHaveNewVertex[c_index]);
            }
        }//for (int i = 0; i < 3; ++i)

    }//for (int face_id = 0; face_id < oldFaceVerticesLength; face_id += 3)
    //calculate all old vertexes' new coordinates.
    std::vector<Cartesian3> updateVertices;
    for (int i = 0; i < oldVerticesNum; ++i) {
        updateVertices.push_back(adjustOldVertexCoor(i));
    }
    //update all old vertex coordinates with new coordinates.
    for (int i = 0; i < oldVerticesNum; ++i) {
        vertices[i] = updateVertices[i];
    }
    //update new center face vertexes to faceVertices.
    for (int i = 0; i < oldFaceVerticesLength; ++i) {
        faceVertices[i] = newCenterFaces[i];
    }
    //the last step, generating the new data of directed edges.
    generateDirectedEdge();
    //the optional step, calculating new normals
    generateNewNormals();
}

/***
 * this function shoulders generating new data of directed edges
 */
void DirectedEdgeSurface::generateDirectedEdge()
{
    //clear the old data, as we need to generate all new first directed edges and other half edges.
    firstDirectedEdge.clear();
    firstDirectedEdge.resize(vertices.size(), NOT_FOUND);
    otherHalf.clear();
    otherHalf.resize(faceVertices.size(), NOT_FOUND);

    int currentFaceVerticesLength = faceVertices.size();

    for(int edgeId = 0;edgeId < currentFaceVerticesLength;edgeId++)
    {
        int faceId = edgeId / 3;
        //from ------> end
        //follow the rule of the provided files' directed edge.
        //the first directed edge follow the end vertex, not from one.
        int endVertex = faceVertices[edgeId];
        int fromVertex = faceVertices[getPreviousEdge(edgeId)];

        if(fromVertex != endVertex)
        {
            if(firstDirectedEdge[fromVertex] == NOT_FOUND){
                //the current edge id is belonging to the the from vertex.
                firstDirectedEdge[fromVertex] = edgeId;
            }

            if(otherHalf[edgeId] == NOT_FOUND){
                //when we look for the other half edge, we just need check the later triangles.
                //as previous have done that and other half won't appear in the current triangle.
                for(int nextCandidateEdgeId = 3 * (faceId + 1); nextCandidateEdgeId < currentFaceVerticesLength; nextCandidateEdgeId++)
                {
                    int nextEndVertex = faceVertices[nextCandidateEdgeId];
                    int nextFromVertex = faceVertices[getPreviousEdge(nextCandidateEdgeId)];

                    if(endVertex == nextFromVertex && fromVertex == nextEndVertex)
                    {
                        otherHalf[edgeId] = nextCandidateEdgeId;
                        otherHalf[nextCandidateEdgeId] = edgeId;
                    }
                }//for(int nextCandidateEdgeId = 3 * (faceId + 1);
            }//if(otherHalf[edgeId] == NOT_FOUND)
        }//if(fromVertex != endVertex)
    }//for(int edgeId = 0;edgeId < currentFaceVerticesLength;edgeId++)
}
/***
 * adjust the old vertex coordinates
 * @param vertexId is vertex index
 * @return adjusted coordinates
 */
Cartesian3 DirectedEdgeSurface::adjustOldVertexCoor(int vertexId)
{
//      c    b
//     /\----/\
//    /  \  /  \
//   /    \/    \
//  d------a-----g
//   \    /\    /
//    \  /  \  /
//     \/____\/
//      e     f
// n is the degree of a
// a = 1/n[5/8 - (3/8 + 1/4 * cos(2pi/n)) ^2]
// another formula is : if n==3, the weight is 3.0/16.0, otherwise is 3.0/(8.0*n);
// so I choose to use the later one, which is simpler. I have checked that, they work equally.

    const static float coefficient_0 = 5.0/8.0;
    const static float coefficient_1 = 3.0/8.0;
    const static float coefficient_2 = 1.0/4.0;

    Cartesian3 sum = {0,0,0};
    int n = 0;

    //find all old directly connected vertexes
    int edgeId = firstDirectedEdge[vertexId];
    do
    {
        //accumulate all vertices' coordinates;
        sum = sum + vertices[faceVertices[getNextEdge(edgeId)]];
        edgeId = getNextEdge(otherHalf[edgeId]);
        n++;
    }
    while (edgeId != firstDirectedEdge[vertexId]);
    float weight = 3.0/16.0;
    if (n != 3)
    {
        weight = 3.0/(8.0*n);
    }
    //calculate the weight of the formula
//  float weight = (1.0 / n * (coefficient_0 - std::pow(coefficient_1 + coefficient_2 * std::cos((2 * M_PI) / n), 2)));
    return vertices[vertexId] * (float)(1 - n * weight) + sum * weight;
}

/***
 * get the face id of the current edge's other half edge
 * @param index, the current edge's index
 * @return a face id
 */
int DirectedEdgeSurface::getFaceId_opposite_theEdge(int index)
{
    return otherHalf[index]/3;
}

int DirectedEdgeSurface::getNextEdge(int index)
{
    //easy one, get the next edge inside a triangle.
    return (index % 3) == 2 ? index - 2 : index + 1;
}

int DirectedEdgeSurface::getPreviousEdge(int index)
{
    //easy one, get the previous edge inside a triangle.
    return (index % 3) != 0 ? index - 1 : index + 2;
}

Cartesian3& DirectedEdgeSurface::getVertexCoor(int32_t index)
{
    //just get the vertex coordinates
    return vertices[faceVertices[index]];
}

/***
 * save the file
 */
void DirectedEdgeSurface::saveCurrentData()
{
    string newFileName = generateNewFileName();
    this->newFileName = newFileName;
    std::ofstream out(newFileName);
    WriteObjectStream(out);
    out.close();
}
/***
 * generate new normals
 */
void DirectedEdgeSurface::generateNewNormals()
{
    normals.resize(vertices.size(),Cartesian3());
    for (int i = 0; i < faceVertices.size(); i+=3)
    {
        int a_index = faceVertices[i];
        int b_index = faceVertices[i+1];
        int c_index = faceVertices[i+2];

        Cartesian3 normal = (vertices[b_index] - vertices[a_index]).cross((vertices[c_index]-vertices[a_index]));
        normals[a_index] =  normals[a_index] = normal;
        normals[b_index] =  normals[b_index] = normal;
        normals[c_index] =  normals[c_index] = normal;
    }

    for (int i = 0; i < normals.size(); ++i) {
        Cartesian3 currentNormal = normals[i];
        float length = sqrt(currentNormal.x * currentNormal.x + currentNormal.y*currentNormal.y + currentNormal.z*currentNormal.z);
        if (length == 0.0)
        {
            length = 1.0;
        }
        normals[i] = currentNormal / length;
    }

}

/*!
 * generate the new file name
 * @return the name of new file.
 * the form of new file name is tetrahedron_sub_division_by_2.diredgenormal,
 * by_2 means the file is divided by 2 , if 4 means by divided by 4 compared to the original file.
 */
std::string DirectedEdgeSurface::generateNewFileName()
{
    string newFileDir = strdup("./resultFiles/");
    string newFileSuffix = strdup(".diredgenormal");
    string newFileName = string();
    string fileFullPath = this->filename;
    string::size_type targetPosition = fileFullPath.find_last_of("/") + 1;
    string file = fileFullPath.substr(targetPosition, fileFullPath.length() - targetPosition);

    string objName = file.substr(0, file.rfind("."));

    if (access(newFileDir.c_str(), 0) == -1)	{
        mkdir(newFileDir.c_str(), 0777);
    }

    newFileName = newFileDir + objName + "_sub_division_by_" + std::to_string(this->subFactor) + newFileSuffix;
    return newFileName;
}