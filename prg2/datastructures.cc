// Datastructures.cc
//
// Student name: 
// Student email: 
// Student number: 

#include "datastructures.hh"

#include <random>

#include <cmath>
std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

// Modify the code below to implement the functionality of the class.
// Also remove comments from the parameter names when you implement
// an operation (Commenting out parameter name prevents compiler from
// warning about unused parameters on operations you haven't yet implemented.)

Datastructures::Datastructures()
{
    allAffiliations = std::vector<Affiliation>();
    allAffiliationIDs = std::vector<AffiliationID>();
    //allAffiliations.reserve(100000);
    affiliationByIds = std::unordered_map<AffiliationID, Affiliation>();
    //allAffiliationName_Ids.reserve(100000);
    //allPublications = std::vector<Publication>();
    //allPublications.reserve(100000);
    publicationByIds = std::unordered_map<PublicationID, Publication>();
    affiliationToPublications = std::unordered_map<AffiliationID, std::vector<PublicationID>>();
}

Datastructures::~Datastructures()
{
    clear_all();
}

unsigned int Datastructures::get_affiliation_count()
{
    // Replace the line below with your implementation
    return allAffiliations.size();
}

void Datastructures::clear_all()
{
    // Replace the line below with your implementation
    allAffiliations.clear();
    allAffiliationIDs.clear();
    affiliationByIds.clear();
    //allPublications.clear();
    publicationByIds.clear();
    affiliationToPublications.clear();

    std::unordered_set<Connection*> connectionSet;

    // Copy pointers from allConnections to the set
    for (Connection* connPtr : allConnections) {
        connectionSet.insert(connPtr);
    }

    // Copy pointers from connectionById to the set
    for (const auto& kvp : connectionById) {
        const std::vector<Connection*>& connections = kvp.second;
        for (Connection* connPtr : connections) {
            connectionSet.insert(connPtr);
        }
    }
    for (Connection* connPtr : connectionSet) {
        delete connPtr; // Deallocate memory for each Connection object
    }
    connectionSet.clear();
    allConnections.clear();
    connectionById.clear();

}

std::vector<AffiliationID> Datastructures::get_all_affiliations()
{
    return allAffiliationIDs;
}

bool Datastructures::add_affiliation(AffiliationID id, const Name &name, Coord xy)
{
    double distance = xy.x * xy.x + xy.y * xy.y;
    Affiliation new_affiliation = {id, name, xy, distance};
    affiliationByIds.insert({id, new_affiliation});
    allAffiliations.push_back(new_affiliation);
    allAffiliationIDs.push_back(id);
    affiliationToPublications[id] = {};
    //affiliationsByName.push_back(std::make_pair(name, id));
    connectionById[id] = std::vector<Connection*>();
    return true;

}

Name Datastructures::get_affiliation_name(AffiliationID id)
{
    if (affiliationByIds.find(id) != affiliationByIds.end()) {
        return affiliationByIds[id].name;
    }
    return NO_NAME;
}

Coord Datastructures::get_affiliation_coord(AffiliationID id)
{
    if (affiliationByIds.find(id) != affiliationByIds.end()) {
        return affiliationByIds[id].xy;
    }
    return NO_COORD;
}

std::vector<AffiliationID> Datastructures::get_affiliations_alphabetically()
{
    std::vector<AffiliationID> result(allAffiliations.size());

    // Sort the affiliations vector alphabetically based on names
    std::sort(allAffiliations.begin(), allAffiliations.end(),
              [](const Affiliation &a, const Affiliation &b) {
                  return a.name < b.name;
              });

    for (std::size_t i = 0; i < allAffiliations.size(); ++i) {
        result[i] = allAffiliations[i].id;
    }
    return result;

}

std::vector<AffiliationID> Datastructures::get_affiliations_distance_increasing()
{
    // Sort the affiliations vector based on coordinates
    std::sort(allAffiliations.begin(), allAffiliations.end(),
              [](const Affiliation &a, const Affiliation &b) {
                  if (a.distance != b.distance) {
                      return a.distance < b.distance;
                  } else {
                      // If distances are equal, compare y coordinates
                      return a.xy.y < b.xy.y;
                  }
              });
    std::vector<AffiliationID> result(allAffiliations.size());
    for (std::size_t i = 0; i < allAffiliations.size(); ++i) {
        result[i] = allAffiliations[i].id;
    }

    return result;
}

AffiliationID Datastructures::find_affiliation_with_coord(Coord xy)
{
    std::sort(allAffiliations.begin(), allAffiliations.end(), [](const Affiliation& a, const Affiliation& b) {
        return a.xy.x < b.xy.x || (a.xy.x == b.xy.x && a.xy.y < b.xy.y);
    });

    // Perform binary search to find the affiliation with the specified Coord
    auto it = std::lower_bound(allAffiliations.begin(), allAffiliations.end(), xy,
                               [](const Affiliation& a, const Coord& xy) {
                                   return a.xy.x < xy.x || (a.xy.x == xy.x && a.xy.y < xy.y);
                               });

    if (it != allAffiliations.end() && it->xy.x == xy.x && it->xy.y == xy.y) {
        return it->id; // Return the AffiliationID if found
    } else {
        return NO_AFFILIATION;
    }
}

bool Datastructures::change_affiliation_coord(AffiliationID id, Coord newcoord)
{
    double distance = newcoord.x * newcoord.x + newcoord.y * newcoord.y;

    auto it = std::find_if(allAffiliations.begin(), allAffiliations.end(),
                           [id](const Affiliation &affiliation) {
                               return affiliation.id == id;
                           });

    if (it != allAffiliations.end()) {
        // If the affiliation is found, update its coordinates, and distance
        it->xy = newcoord;
        it->distance = distance;
    }
    else {
        return false;
    }

    if (affiliationByIds.find(id) != affiliationByIds.end()) {
        affiliationByIds[id].xy = newcoord;
        affiliationByIds[id].distance = distance;
        return true;
    }
    return false;
}

bool Datastructures::add_publication(PublicationID id, const Name &name, Year year, const std::vector<AffiliationID> &affiliations)
{
    Publication new_publication = {id, name, year, affiliations};
    publicationByIds.emplace(id, new_publication);
    for (const auto& affid:affiliations) {
        if (affiliationToPublications.find(affid) == affiliationToPublications.end()) {
            affiliationToPublications.insert({affid, std::vector<PublicationID>()});
        }
        affiliationToPublications[affid].push_back(id);
    }
    create_connection(new_publication);
    return true;
}

std::vector<PublicationID> Datastructures::all_publications()
{
    std::vector<PublicationID> result;
    result.reserve(publicationByIds.size());

    for (auto it = publicationByIds.begin(); it != publicationByIds.end(); ++it) {
        result.push_back(it->first);
    }
    return result;

}

Name Datastructures::get_publication_name(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        return publicationByIds[id].name;
    }
    return NO_NAME;

}

Year Datastructures::get_publication_year(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        return publicationByIds[id].year;
    }
    return NO_YEAR;
}

std::vector<AffiliationID> Datastructures::get_affiliations(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        return publicationByIds[id].by_affiliations;
    }
    return {NO_AFFILIATION};
}

bool Datastructures::add_reference(PublicationID id, PublicationID parentid)
{
    auto referencing = publicationByIds.find(parentid);
    auto referenced = publicationByIds.find(id);
    if (referencing != publicationByIds.end() && referenced != publicationByIds.end()) {
        referenced->second.parent = parentid;
        referencing->second.children.push_back(id);
        return true;
    }
    return false;
}

std::vector<PublicationID> Datastructures::get_direct_references(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        return publicationByIds[id].children;
    }
    return {};
}

bool Datastructures::add_affiliation_to_publication(AffiliationID affiliationid, PublicationID publicationid)
{
    auto pub = publicationByIds.find(publicationid);
    if (pub != publicationByIds.end() && affiliationByIds.find(affiliationid) != affiliationByIds.end()) {
        pub->second.by_affiliations.push_back(affiliationid);
        affiliationToPublications[affiliationid].push_back(publicationid);
        create_connection(pub->second, affiliationid);
        return true;
    }
    return false;

}

std::vector<PublicationID> Datastructures::get_publications(AffiliationID id)
{
    if (affiliationToPublications.find(id) != affiliationToPublications.end()) {
        return affiliationToPublications[id];
    }
    return {NO_PUBLICATION};

}

PublicationID Datastructures::get_parent(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        return publicationByIds[id].parent;
    }
    return NO_PUBLICATION;
}

std::vector<std::pair<Year, PublicationID> > Datastructures::get_publications_after(AffiliationID affiliationid, Year year)
{

    if (affiliationToPublications.find(affiliationid) != affiliationToPublications.end()) {
        std::vector<std::pair<Year, PublicationID>> result;
        auto aff = affiliationToPublications[affiliationid];
        for (auto pubid : aff) {
            if (publicationByIds.find(pubid) != publicationByIds.end()) {
                if (publicationByIds[pubid].year >= year) {
                    result.emplace_back(publicationByIds[pubid].year, pubid);
                }
            }
        }
        // Sort the publications in ascending order by publication year and, if necessary, by ID
        std::sort(result.begin(), result.end(), [](const auto &a, const auto &b) {
            if (a.first == b.first) {
                return a.second < b.second;
            }
            return a.first < b.first;
        });

        return result;
    }
    return {{NO_YEAR, NO_PUBLICATION}};
}

std::vector<PublicationID> Datastructures::get_referenced_by_chain(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        auto current = publicationByIds[id];
        std::vector<PublicationID> result;
        while (current.parent != NO_PUBLICATION) {
            result.push_back(current.parent);
            current = publicationByIds[current.parent];
        }
        return result;
    }
    return {NO_PUBLICATION};
}


std::vector<PublicationID> Datastructures::get_all_references(PublicationID id)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        auto pub = publicationByIds[id];
        std::vector<PublicationID> result;
        pretraversal(id, result);
        if (result.size() > 1) {
            result.erase(result.begin());
            return result;
        }
        return {};
    }
    else {
        return {NO_PUBLICATION};
    }
}

void Datastructures::pretraversal(PublicationID id, std::vector<PublicationID>&result)
{
    if (publicationByIds.find(id) != publicationByIds.end()) {
        result.push_back(id);
        for (PublicationID child : publicationByIds[id].children) {
            pretraversal(child, result);
        }
    }
}

std::vector<AffiliationID> Datastructures::get_affiliations_closest_to(Coord xy)
{
    // Create a vector of pairs (AffiliationID, Affiliation)
    std::vector<std::pair<AffiliationID, Affiliation>> sortedPairs(affiliationByIds.begin(), affiliationByIds.end());

    // Sort the vector based on the values (Affiliation)
    std::sort(sortedPairs.begin(), sortedPairs.end(), [xy](const auto& a, const auto& b) {

        double distanceA = std::sqrt(std::pow(xy.x - a.second.xy.x, 2) + std::pow(xy.y - a.second.xy.y, 2));
        double distanceB = std::sqrt(std::pow(xy.x - b.second.xy.x, 2) + std::pow(xy.y - b.second.xy.y, 2));
        if (distanceA == distanceB) {
            return a.second.xy.y < b.second.xy.y;
        }
        return distanceA < distanceB;

    });

    std::vector<AffiliationID> result;
    result.reserve(3);
    for (const auto& affiliation : sortedPairs) {
        if (result.size() >= 3) {
            break;
        }
        result.push_back(affiliation.second.id);
    }
    return result;

}

bool Datastructures::remove_affiliation(AffiliationID id)
{
    // Delete affiliation from allAffiliations
    auto affit = std::find_if(allAffiliations.begin(), allAffiliations.end(),
                              [id](const Affiliation &affiliation) {
                                  return affiliation.id == id;
                              });

    if (affit != allAffiliations.end()) {
        allAffiliations.erase(affit);
    }

    // Delete affiliation ID ffrom allAffiliationIDs
    auto idt = find(allAffiliationIDs.begin(), allAffiliationIDs.end(),
                    id);

    if (idt != allAffiliationIDs.end()) {
        allAffiliationIDs.erase(idt);
    }
    // Delete affiliation from affiliationByIds
    auto it = affiliationByIds.find(id);
    if (it != affiliationByIds.end()) {
        // Delete affiliation id from affiliationToPublications and delete the
        // record of affiliation id from each of its publications
        auto it2 = affiliationToPublications.find(id);
        if (it2 != affiliationToPublications.end()) {
            for (auto pubId : it2->second) {
                auto pubIt = publicationByIds.find(pubId);
                if (pubIt != publicationByIds.end()) {
                    // Assuming by_affiliations is a std::vector<AffiliationID>
                    pubIt->second.by_affiliations.erase(std::remove(pubIt->second.by_affiliations.begin(),
                                                                    pubIt->second.by_affiliations.end(),
                                                                    id),
                                                        pubIt->second.by_affiliations.end());
                }
            }
            affiliationToPublications.erase(it2);
        }
        affiliationByIds.erase(it);
        return true;
    }
    return false;

}

PublicationID Datastructures::get_closest_common_parent(PublicationID id1, PublicationID id2)
{
    std::vector<PublicationID> parents1 = get_referenced_by_chain(id1);
    std::vector<PublicationID> parents2 = get_referenced_by_chain(id2);
    for (const PublicationID& id : parents2) {
        auto publicationIt = std::find_if(parents1.begin(), parents1.end(),
                                          [id](const PublicationID& publicationId) {
                                              return publicationId == id;
                                          });
        if (publicationIt != parents1.end()) {
            // Found the first common element
            return id;
        }
    }
    return NO_PUBLICATION;
}

bool Datastructures::remove_publication(PublicationID publicationid)
{
    auto it = publicationByIds.find(publicationid);
    if (it != publicationByIds.end()) {
        // Remove references to this publication from the parent publication
        if (it->second.parent != NO_PUBLICATION) {
            auto parentIt = publicationByIds.find(it->second.parent);

            auto removeIt = std::remove(parentIt->second.children.begin(), parentIt->second.children.end(), publicationid);
            parentIt->second.children.erase(removeIt, parentIt->second.children.end());
        }

        // Delete the parent id from this publication's children
        if (!it->second.children.empty()) {
            for (PublicationID child : it->second.children) {
                auto childIt = publicationByIds.find(child);
                childIt->second.parent = NO_PUBLICATION;
            }
        }

        if (!it->second.by_affiliations.empty()) {
            for (AffiliationID affiliation : it->second.by_affiliations) {
                auto affiliationIt = affiliationToPublications.find(affiliation);
                if (affiliationIt != affiliationToPublications.end()) {
                    auto remove_it = std::find(affiliationIt->second.begin(), affiliationIt->second.end(), publicationid);
                    if (remove_it != affiliationIt->second.end()) {
                        affiliationIt->second.erase(remove_it);
                    }
                }
            }
        }
        // Remove the publication from the list of all publications
        publicationByIds.erase(publicationid);
        return true;
    }
    return false;

}

std::vector<PublicationID> Datastructures::get_common_publication(AffiliationID id1, AffiliationID id2)
{
    std::vector<PublicationID> pub1 = get_publications(id1);
    std::vector<PublicationID> pub2 = get_publications(id2);

    std::vector<PublicationID> commonPublications;
    if (id1 == id2) {
        return commonPublications;
    }
    for (const PublicationID id : pub1) {
        if (std::find(pub2.begin(), pub2.end(), id) != pub2.end()) {
            commonPublications.push_back(id);
        }
    }

    return commonPublications;
}
/*
void Datastructures::create_connection_for_all()
{
    for (const auto& id : allAffiliationIDs) {
        create_connection(id);
    }

}
*/
void Datastructures::create_connection(Publication pub, AffiliationID aff_to_fix)
{
    if (aff_to_fix == NO_AFFILIATION) {
        for (size_t i = 0; i < pub.by_affiliations.size(); ++i) {

            for (size_t j = i + 1; j < pub.by_affiliations.size(); ++j) {
                AffiliationID target;
                AffiliationID source;
                if (pub.by_affiliations[i] < pub.by_affiliations[j]) {
                    source = pub.by_affiliations[i];
                    target = pub.by_affiliations[j];
                }
                else {
                    source = pub.by_affiliations[j];
                    target = pub.by_affiliations[i];
                }
                if (target == source) {
                    break;
                }
                bool exist = false;
                auto it = connectionById.find(source);

                if (it != connectionById.end()) {
                    std::vector<Connection*>& connections = it->second;
                    for (Connection* conn : connections) {
                        if (conn->aff1 == source && conn->aff2 == target) {
                            // Found a connection with the same source and target
                            exist = true;
                            conn->weight += 1;
                            break;
                        }
                    }

                    if (!exist) {
                        Connection* connection = new Connection{source, target, 1};
                        Connection* connection_reverse = new Connection{target, source, 1};
                        allConnections.push_back(connection);

                        connectionById[source].push_back(connection);
                        connectionById[target].push_back(connection_reverse);

                    }

                }
                if (exist) {
                    auto it_reverse = connectionById.find(target);
                    if (it_reverse != connectionById.end()) {
                        std::vector<Connection*>& connections_reverse = it_reverse->second;
                        for (Connection* conn : connections_reverse) {
                            if (conn->aff1 == target && conn->aff2 == source) {
                                conn->weight += 1;
                                break;
                            }
                        }
                    }

                }
            }
        }
    } else {
        for (size_t i = 0; i < pub.by_affiliations.size(); ++i) {
            AffiliationID target;
            AffiliationID source;
            if (pub.by_affiliations[i] < aff_to_fix) {
                source = pub.by_affiliations[i];
                target = aff_to_fix;
            }
            else {
                source = aff_to_fix;
                target = pub.by_affiliations[i];
            }
            if (target == source) {
                break;
            }
            bool exist = false;
            auto it = connectionById.find(source);

            if (it != connectionById.end()) {
                std::vector<Connection*>& connections = it->second;
                for (Connection* conn : connections) {
                    if (conn->aff1 == source && conn->aff2 == target) {
                        // Found a connection with the same source and target
                        exist = true;
                        conn->weight += 1;
                        break;
                    }
                }

                if (!exist) {
                    Connection* connection = new Connection{source, target, 1};
                    Connection* connection_reverse = new Connection{target, source, 1};
                    allConnections.push_back(connection);

                    connectionById[source].push_back(connection);
                    connectionById[target].push_back(connection_reverse);

                }

            }
            if (exist) {
                auto it_reverse = connectionById.find(target);
                if (it_reverse != connectionById.end()) {
                    std::vector<Connection*>& connections_reverse = it_reverse->second;
                    for (Connection* conn : connections_reverse) {
                        if (conn->aff1 == target && conn->aff2 == source) {
                            conn->weight += 1;
                            break;
                        }
                    }
                }

            }
        }
    }

    /*
    auto source = affiliationByIds.find(id);
    if (source != affiliationByIds.end()) {

        for (const auto& target : allAffiliations) {
            std::vector<PublicationID> commonPub = get_common_publication(id, target.id);
            int weight = int(commonPub.size());
            if (weight != 0) {
                // double length = std::sqrt(std::pow(source->second.xy.x -target.xy.x, 2) + std::pow(source->second.xy.y - target.xy.y, 2));
                Connection conn = {id, target.id, weight};


                bool add_to_all_connections = true;
                for (auto&c : allConnections) {
                    if ((c.aff1 == conn.aff1 && c.aff2 == conn.aff2) || (c.aff1 == conn.aff2 && c.aff2 == conn.aff1)) {
                        add_to_all_connections = false;
                        if (c.weight != weight) {
                            c.weight = weight;
                        }
                    }
                }
                if (add_to_all_connections) {
                    if (id > target.id) {
                        Connection conn2 = {target.id, id, weight};
                        allConnections.push_back(conn2);
                    } else {
                        allConnections.push_back(conn);
                    }

                }

                if (connectionById.find(id) == connectionById.end()) {
                    std::vector<Connection> newVector;
                    connectionById.insert({id, newVector});
                }
                bool connectionExists = false;
                // bool updateWeight = false;
                for (auto& connection : connectionById[id]) {
                    if (connection == conn) {
                        connectionExists = true;
                        break;
                    }
                    if (connection.aff1 == conn.aff1 && connection.aff2 == conn.aff2) {
                        connectionExists = true;
                        connection.weight = weight;
                        break;
                    }
                }

                // If the connection doesn't exist for this 'id', push it back
                if (!connectionExists) {
                    connectionById[id].push_back(conn);
                }

            }
        }

    }
    */
}

std::vector<Connection> Datastructures::get_connected_affiliations(AffiliationID id)
{
    //create_connection_for_all();
    //create_connection(id);

    // std::vector<Connection> result;
    if (connectionById.find(id) != connectionById.end()) {
        std::vector<Connection> result;
        result.reserve(connectionById[id].size());
        for (const auto& ptr : connectionById[id]) {
            if (ptr != nullptr) {
                result.emplace_back(*ptr); // Dereference pointer and push_back into allConnections
            }
        }
        return result;
    }

    return {};
}

std::vector<Connection> Datastructures::get_all_connections()
{
    // create_connection_for_all();
    std::vector<Connection> result;
    result.reserve(allConnections.size());

    // Convert pointers to objects
    for (const auto& ptr : allConnections) {
        if (ptr != nullptr) {
            result.emplace_back(*ptr); // Dereference pointer and push_back into allConnections
        }
    }
    return result;
    //return allConnections;
}

bool Datastructures::findPath(AffiliationID current, AffiliationID target, Path& path, std::unordered_map<AffiliationID, bool>& visited) {

    if (current == target) {
        return true;
    }

    visited[current] = true;

    for (const Connection* connection : connectionById[current]) {

        // Check if the connection is not visited and if it leads to the target
        if (!visited[connection->aff2]) {
            path.push_back(*connection);

            if (findPath(connection->aff2, target, path, visited)) {

                return true;
            }
            path.pop_back();
        }
    }

    return false;
}

Path Datastructures::get_any_path(AffiliationID source, AffiliationID target)
{

    // create_connection_for_all();
    Path path;
    std::unordered_map<AffiliationID, bool> visited;

    if (source == target) {
        return path; // Empty path as source and target are the same
    }

    for (Connection* conn : allConnections) {
        visited[conn->aff1] = false;
        visited[conn->aff2] = false;
    }

    if (findPath(source, target, path, visited)) {
        return path;
    }

    return {}; // Empty path if no path found


}

Path Datastructures::get_path_with_least_affiliations(AffiliationID source, AffiliationID target)
{

    // create_connection_for_all();
    std::unordered_set<AffiliationID> visited;
    std::queue<Path> paths;
    paths.push({});

    while (!paths.empty()) {
        std::vector<Connection> currentPath = paths.front();
        paths.pop();

        Connection lastConnection;
        if (!currentPath.empty()) {
            lastConnection = currentPath.back();
        }

        AffiliationID currentAffiliation = (currentPath.empty()) ? source : lastConnection.aff2;

        if (currentAffiliation == target) {
            return currentPath;
        }

        visited.insert(currentAffiliation);

        for (const auto& conn : connectionById[currentAffiliation]) {
            if (visited.find(conn->aff2) == visited.end()) {
                std::vector<Connection> newPath = currentPath;
                newPath.push_back(*conn);
                paths.push(newPath);
            }
        }
    }

    return {}; // Empty path if no path found
}

double Datastructures::calculateFriction(const Path& path) {
    double minWeight = std::numeric_limits<double>::max();
    for (const auto& conn : path) {
        if (conn.weight < minWeight) {
            minWeight = conn.weight;
        }
    }

    return 1.0 / minWeight;
}

void Datastructures::findAllPaths(AffiliationID source, AffiliationID target, Path& currentPath, std::vector<Path>& allPaths, std::unordered_set<AffiliationID>& visited) {

    visited.insert(source);

    if (source == target) {
        // If we reach the target, store the current path in allPaths
        allPaths.push_back(currentPath);
    } else {
        for (const auto& conn : connectionById[source]) {
            if (visited.find(conn->aff2) == visited.end()) {
                currentPath.push_back(*conn);
                findAllPaths(conn->aff2, target, currentPath, allPaths, visited);
                currentPath.pop_back();
            }
        }
    }

    visited.erase(source);

}

Path Datastructures::get_path_of_least_friction(AffiliationID source, AffiliationID target)
{
    // create_connection_for_all();
    Path path;
    std::vector<Path> allPaths = std::vector<Path>();
    std::unordered_set<AffiliationID> visited;
    findAllPaths(source, target, path, allPaths, visited);
    double minFriction = std::numeric_limits<double>::max();
    int minPathLength = std::numeric_limits<int>::max();
    Path bestPath;
    for (auto& path : allPaths) {
        double friction = calculateFriction(path);
        if (friction < minFriction) {
            minFriction = friction;
            minPathLength = path.size();
            bestPath = path;
        }
        if (friction == minFriction) {
            if (int(path.size()) < minPathLength) {
                minPathLength = path.size();
                bestPath = path;
            }
        }
    }
    return bestPath;
}

double Datastructures::calculateConnectionLength(const Connection conn) {
    auto source = affiliationByIds.find(conn.aff1);
    auto target = affiliationByIds.find(conn.aff2);
    double length = 0;
    if (source != affiliationByIds.end() && source != affiliationByIds.end()) {
        length = std::sqrt(std::pow(source->second.xy.x -target->second.xy.x, 2) + std::pow(source->second.xy.y - target->second.xy.y, 2));

    }
    return length;
}

std::pair<double, PathWithDist> Datastructures::Dijkstra_shortest(AffiliationID source, AffiliationID target) {
    std::priority_queue<std::pair<int, AffiliationID>, std::vector<std::pair<int, AffiliationID>>, std::greater<std::pair<int, AffiliationID>>> pq;
    std::unordered_map<AffiliationID, int> distance;
    std::unordered_map<AffiliationID, Connection*> previous; // Store previous connections

    // Initialize distances to infinity
    for (auto& kp : connectionById) {
        distance[kp.first] = std::numeric_limits<int>::max();
    }

    distance[source] = 0;
    pq.push({0, source});

    while (!pq.empty()) {
        AffiliationID current = pq.top().second;
        int currentDistance = pq.top().first;
        pq.pop();

        if (current == target) {
            PathWithDist shortestPath;
            while (previous.find(current) != previous.end()) {
                shortestPath.push_back(std::make_pair(*previous[current], calculateConnectionLength(*previous[current])));
                current = previous[current]->aff1; // Move to the previous node in the path
            }
            std::reverse(shortestPath.begin(), shortestPath.end()); // Reverse to get correct order
            return { distance[target], shortestPath }; // Return shortest path length and the path
        }

        for (Connection* conn : connectionById[current]) {

            double connection_length = calculateConnectionLength(*conn);
            int newDistance = currentDistance + connection_length;
            if (newDistance < distance[conn->aff2]) {
                distance[conn->aff2] = newDistance;
                pq.push({newDistance, conn->aff2});
                previous[conn->aff2] = conn; // Record the connection for the shortest path
            }
        }
    }

    return { -1, {} }; // No path found, return -1 for distance and an empty vector for the path
}

PathWithDist Datastructures::get_shortest_path(AffiliationID source, AffiliationID target)
{
    /*
    // create_connection_for_all();
    Path path;
    std::vector<Path> allPaths = std::vector<Path>();
    std::unordered_set<AffiliationID> visited;
    findAllPaths(source, target, path, allPaths, visited);
    double minLength = std::numeric_limits<double>::max();

    PathWithDist bestPath;
    for (auto& path : allPaths) {
        PathWithDist newPath;
        double length = calculateLength(path, newPath);
        if (length < minLength) {
            minLength = length;
            bestPath = newPath;
        }

    }
    return bestPath;
    */
    PathWithDist bestPath;
    double shortestPathLength;


    std::tie(shortestPathLength, bestPath) = Dijkstra_shortest(source, target);
    if (shortestPathLength == -1 || shortestPathLength == 0) {
        return bestPath;
    } else {
        return bestPath;
    }

}



