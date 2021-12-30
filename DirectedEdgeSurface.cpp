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
        geometryStream << "Face " << face << " ";
        
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

// 提纲
// 1,首先 记录 现在的 vertex 数组 size, 因为新产生的定点要直接追加到当前定点数组后面, 所以要知道老 定点究竟是那些
//   对于 normal  firstdirectededge face otherhalf 也是这样记录
// 2, 按 faces 数组进行遍历, 开始生成新的定点,每个边一个新的定点,
//    但是会有的边会重复计算新的定点, 这里用一个新的数据结构去记录一下已经产生过新的定点的边,后面遇到已经标记过的,就跳过去不用重新计算,
//    但是 这个定点要被的索引要被追加到new face 索引表中
// 3, 新生成的定点,直接追加到定点数组后面, 而这些新加入的定点,
//    此时有三个变化了的数组, face, vertex, new face
// 4, 然后按 face中old部分 调整每个老定点的坐标, 并且更新 vertex, 还有face 数组中老 face的定点索引
// 5, 重新计算otherhalf, first directed edge, normal

//          c
//       /    \
//     /        \
//    a - - m - - b
//     \        /
//      \      /
//        \  /
//         d



void DirectedEdgeSurface::loopSubDivision()
{
    int oldVerticesNum = vertices.size();
    int oldFacesNum = faceVertices.size() / 3;

    int oldFaceVerticesLength = faceVertices.size();
    std::vector<int > newCenterFaces;
    std::vector<int> edgesHaveNewVertex;
    edgesHaveNewVertex.resize(otherHalf.size(),NOT_FOUND);

    for (int face_id = 0; face_id < oldFacesNum; face_id++) {
        for (int i = 0; i < 3; ++i) {
            //inside a face, 这里面三个索引也是对应边的索引
            int a_index = face_id * 3 + i % 3;
            int b_index = face_id * 3 + (i + 1) % 3;
            int c_index = face_id * 3 + (i + 2) % 3;

            int b_other_half_edge = otherHalf[b_index];
            if (edgesHaveNewVertex[b_index] == NOT_FOUND && edgesHaveNewVertex[b_other_half_edge] == NOT_FOUND)
            {
                //如果这条边,没有生成过新的定点,那么生成新的定点

                int d_index = getNextEdge(b_other_half_edge);
//            int d_index = getFaceId_opposite_theEdge(b_index) * 3 + (getNextEdge(b_other_half_edge) % 3);

                Cartesian3 vertex_d_coor;
                //vertex d 的坐标
                vertex_d_coor = 0.375f * (getVertexCoor(a_index) + getVertexCoor(b_index))
                                + 0.125 * (getVertexCoor(c_index) + getVertexCoor(d_index));

                vertices.push_back(vertex_d_coor);
                //一对边都得存一下 插入点的索引,她们公用同一个点
                int newVertexIndex = vertices.size() - 1;
                edgesHaveNewVertex[b_index] = newVertexIndex;
                edgesHaveNewVertex[b_other_half_edge] = newVertexIndex;
                newCenterFaces.push_back(newVertexIndex);
                faceVertices.push_back(newVertexIndex);

            }//if (edgesHaveNewVertex[b_index] != NOT_FOUND && edgesHaveNewVertex[b_other_half_edge] != NOT_FOUND)
            else {
                //如果这条边已经生成新的定点了,
                int newVertexIndex = edgesHaveNewVertex[b_index];
                newCenterFaces.push_back(newVertexIndex);
                faceVertices.push_back(newVertexIndex);
            }
            if (faceVertices.size() % 3 == 0)
            {
                faceVertices.push_back(faceVertices[faceVertices.size()-1]);
            }
            faceVertices.push_back(b_index);
            if (i==2){
                //当前三角形的最后一次,要补充最后一个三角形的定点,在构造新的三角形的时候,不然会少一个顶点
                faceVertices.push_back(edgesHaveNewVertex[c_index]);
            }
        //到这里,就搞定了每个三角形中插入新的顶点,并且构造了新的face,全部新顶点的一个,另外三个,
        }//for (int i = 0; i < 3; ++i)



    }//for (int face_id = 0; face_id < oldFaceVerticesLength; face_id += 3)
    //接下来,更新所有的 老顶点的坐标
    for (int i = 0; i < oldVerticesNum; ++i) {
        vertices[i] = adjustOldVertexCoor(i);
    }
    //更新old face 的顶点索引
    for (int i = 0; i < oldFaceVerticesLength; ++i) {
        faceVertices[i] = newCenterFaces[i];
    }
    //最后一步,生成 directed edge
    generateDirectedEdge();
}


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
        //遵循和 当前文件一样的规则, 以 end 顶点索引为 edge索引
        int endVertex = faceVertices[edgeId];
        int fromVertex = faceVertices[getPreviousEdge(edgeId)];

        if(fromVertex != endVertex)
        {
            if(firstDirectedEdge[endVertex] == NOT_FOUND){
                firstDirectedEdge[endVertex] = getPreviousEdge(edgeId);
            }

            if(otherHalf[edgeId] == NOT_FOUND){
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

    const static float coefficient_0 = 5.0/8.0;
    const static float coefficient_1 = 3.0/8.0;
    const static float coefficient_2 = 1.0/4.0;

    Cartesian3 sum = {0,0,0};
    int32_t n = 0;

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

    //calculate the weight of the formula
    float weight = (1.0 / n * (coefficient_0 - std::pow(coefficient_1 + coefficient_2 * std::cos((2 * M_PI) / n), 2)));
    return vertices[vertexId] * (float)(1 - n * weight) + sum * weight;
}

//index : the edge id
int DirectedEdgeSurface::getFaceId_opposite_theEdge(int index)
{
    //当前edge 的 other half edge的 所在的face id
    return otherHalf[index]/3;
}

int DirectedEdgeSurface::getNextEdge(int index)
{
    //每个三角形内,第三条边的id是大于第一条的,所以减去2,此时构成一个环,又回去了
    return (index % 3) == 2 ? index - 2 : index + 1;
}

int DirectedEdgeSurface::getPreviousEdge(int index)
{
    //和next edge一个道理
    return (index % 3) != 0 ? index - 1 : index + 2;
}

Cartesian3& DirectedEdgeSurface::getVertexCoor(int32_t index)
{
    return vertices[faceVertices[index]];
}

//save file
void DirectedEdgeSurface::saveCurrentData()
{
    string newFileName = generateNewFileName();
    this->newFileName = newFileName;
    std::ofstream out(newFileName);
    WriteObjectStream(out);
    out.close();
}


/*!
 * generate the new file name for texture.ppm and normalMap.ppm file
 * @return ResultFileNames, a struct
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