#include<bits/stdc++.h>
#include<unistd.h>
using namespace std;

//maps to store the transactions and blocks created in the network.
map <int,struct Transaction> global_transactions;
map<int,struct Block> global_blocks;
int n;
//block id's start from 1
int global_block_id=1;

int coin_base_transaction_id=0;

//transaction id's start from 101, 0-100 are reserved for the coinbase transactions for 100 blocks which we going to add in blockchain
int txn_id=101;
vector<int> transactions_gen_by_each_node;
vector<int> blocks_gen_by_each_node;
vector<int> blocks_added_by_each_node;
//structure for each transaction
struct Transaction
{
	int id;
	int source;
	int destination;
	double transfer_amount;
	string message;
};
//structure for each block
struct Block
{
	int id;
	vector<Block> pointed_by;  //block which are pointing to this block i.e child's of this block
	int parent_id;
	int created_by;
	double time_created;
	double time_added;
	int length_of_chain;   //length of the chain till this block in which it is added
	
	vector<int> balance_nodes_at_this_block;   //vector to track the balances of the nodes when this block has added to the blockchain
	vector< int> transactions;         // vector to store the transactions present in the block.
	vector<int> transactions_till_now;   // vector to keep track of all the transactions that are there in the chain till this block.
};

//structure to keep track the nodes.
struct Node
{
	int id;
	string type_of_cpu;  
	string type_of_node;
	double amount;
	double hashing_power;
	map < int ,double> transactions;   //map to track the transactions id's seen by this node and the time it has seen that transaction.
	map<int,double> block_ids_seen;    // map to track the block id's seen by this node and the time it has seen that block.
	vector <pair<int,double>> Peers;  //vector to track the connections of this node.
	Block genesis_block;    //every node will has genesis block which is the start of blockchain.
	int working_on;      // variable to keep track the id of block on which the node is currently mining.
};
Node * network;

//struct for the event which will have operation of this event, the id, and the node to which this event is targeted to, and time this event should be handled.
struct Event{
    string operation;
    int id;
    int src;
    double time;
    Event(string operation, int id,int src,double time)
        : operation(operation), id(id),src(src),time(time)
    {
    }
};

//compare_time helps to create the priority_queue having the least time event at the top every time.
struct compare_time {
    bool operator()(Event const& p1, Event const& p2)
    {
        return p1.time > p2.time;
    }
};

//priority_queue having all the events information which should be handled.
priority_queue< Event, vector<Event> ,compare_time> pq;


//find the hashing power of the node which is having slow_cpu based on the no.of high_cpu nodes and no.of low_cpu nodes and hashing power of high_cpu node is 10 times of the low_cpu
double find_hashing_power(double high,double low)
{
	double z=1/( 10 * (high) + low);
	return z;
}

//function to get some random number from an exponential distribution given the mean.
double get_exponential_dist_numb(double mean)
{
	random_device rd;
    mt19937 gen(rd());
    exponential_distribution<> d(1.0/mean);
   	return d(gen);
 
}

//function to add the transaction receive actions for the peers of a node once a node receive a transaction.
void broadcast_transaction(int txn_id,int node, double time)
{
	double cij;
	double dij,latency;
	double received_time=time;
	Transaction txn=global_transactions[txn_id];    //getting the transaction based on the trasanction id.
	network[node].transactions[txn_id]=received_time;  //adding the transaction and the time node has seen it into the transaction map of the node.

	//for every peer of this node.
	for(auto peer: network[node].Peers)
	{
		
		//peer id is peer.first
		// the speed of light for this link is peer.second
		//if the peer has already seen this transaction from other's then don't send this transaction to it.
		if(network[peer.first].transactions.find(txn_id)!=network[peer.first].transactions.end())
		{
		continue;
		}

		//if the peer has not seen this transaction yet calculate the latencies and push the tranasaction receive event targeting the peer.
			if(network[node].type_of_node=="fast" && network[peer.first].type_of_node=="fast")   //if both the node are fast the link speed cij=100Mbps or else cij=5Mbps
			{
				cij=100;
			}
			else
			{
				cij=5;
			}
			//calucating the latencies as mentioned .
			dij=get_exponential_dist_numb( (double) 96.0 / (cij*1000) );

			//latency is size of message by cij. Here we are sending a transaction. It is given size of transaction is 1KB
			latency= ( (double) 8) / (cij*1000);
			//pushing the event transaction receive with id of transaction to sent targetting the peer and the time at which the peer will receive this transaction. 
			//here it.second is Pij speed of light which is randomly choosen from the uniform distribution [10ms to 500 ms] i.e [0.01 to 0.5 secs] for every link
			pq.push(Event("txn_recv",txn_id, peer.first ,time+peer.second+latency+dij));
	}
}

//function to create a transaction and call the broadcast_transaction function to push the transaction_receive actions for the peer of the node which has created this transaction.
void generate_transaction(double time)
{
	
		int source=rand()%n;   //randonly select any node in the network.
		int destination=source;  //randomly select the destination which is not same as source.
		while(source==destination)
		{
			destination=rand()%n;
		}
		//randomly select the amount which needs to be transfered in this transaction.
		int amount=rand()%100;
		//creating the transaction with the choosen values.
		struct Transaction txn;
		txn.source=source;
		txn.destination=destination;
		txn.transfer_amount=amount;
		transactions_gen_by_each_node[source]++;
		string message= "TXN "+to_string(txn_id)+": User: "+to_string(source)+" PAYS "+to_string(amount)+" to User: "+to_string(destination);
		txn.message=message;
		//adding this transaction to global transactino map.
		global_transactions[txn_id]=txn;
		//calling the broadcast_transaction function to push the transaction_recv targetting the peers of this transaction source.
		broadcast_transaction(txn_id,source, time);
		txn_id++;
}

//function to broadcast a block same as broadcast_transaction function.
void broadcast_block(int blk_id,int node, double received_time)
{
	double cij;
	double dij,latency;
	Block blk=global_blocks[blk_id];   //get the block using the block id from the global_blocks.
	network[node].block_ids_seen[blk_id]=received_time;   //add this block the block_ids_seen map of this node and the time at which this node has seen this block.
	//for every peer of this node.
	for(auto it: network[node].Peers)
	{
		//if the peer has not seen this block yet calculate the latencies and push the block_recv actions targetting the peers.
		//the latencies will be calculated same as in broadcast_transaction function.
		if(network[it.first].block_ids_seen.find(blk_id)==network[it.first].block_ids_seen.end())
		{
			if(network[node].type_of_node=="fast" && network[it.first].type_of_node=="fast")
			{
				cij=100;
			}
			else
			{
				cij=5;
			}
			dij=get_exponential_dist_numb( (double) 96.0 / (cij*1000000) );
			//size of the empty block is 1KB
			//size of the non-empty block is 1KB + no.of transactions in the block * size of each transaction i.e 1KB
			latency= ( (double) blk.transactions.size() * 8 + 8 ) / (cij*1000000);
			//push the block_recv action targetting to the peer having the id of block to received, id of peer targetted to, and the time at which the peer will see this block.
			pq.push(Event("block_recv",blk_id,it.first,received_time+it.second+latency+dij));
		}
	}
}
//function to traverse through the graph using dfs.
void dfs(int node,vector<vector<pair<int,double>>> &graph, vector<int> & visited)
{
	for(auto i:graph[node])
	{

		if(visited[i.first]==0)
		{
			visited[i.first]=1;
			dfs(i.first,graph,visited);
		}
	}
}
//function to intialize the network
void initialization( int n, int z0, int z1)
{
	//getting the number of low cpu nodes
	int no_of_low_cpu_type= (int) ( z1 * n/100);
	//finding the hashing power of the node with low cpu type.
	double low_hash_power=find_hashing_power(n- no_of_low_cpu_type, no_of_low_cpu_type );
	//intialize the network
	for(int i=0;i<n;i++)
	{
		network[i].id=i;
		network[i].type_of_cpu="high CPU";
		network[i].type_of_node="fast";
		network[i].amount=100;
		network[i].hashing_power=10*low_hash_power;
		network[i].transactions={};
		network[i].Peers={};
		network[i].working_on=0;
	}
	//getting the no.of nodes which are slow in the network
	int no_of_slow_type= (int) (z0*n/100);
	

	int i=0,random_number;
	//assigning slow for no_of_slow_type nodes randomly
	while(i<no_of_slow_type)
	{
		random_number=rand()%n;
		if(network[random_number].type_of_node=="fast")
		{
			network[random_number].type_of_node="slow";
			i++;
		}
	}
	//assigning low cpy for no_of_low_cpu_type nodes randomly
	i=0;
	while(i<no_of_low_cpu_type)
	{
		random_number=rand()%n;
		if(network[random_number].type_of_cpu=="high CPU")
		{
			network[random_number].type_of_cpu="low CPU";
			network[random_number].hashing_power=low_hash_power;
			i++;
		}
	}

}

//function to create a random graph and check whether it is connected, if not recreate the graph from start.
void graph_creation(int n) {
	
    vector<int> node_ids(n);
    for (int i = 0; i < n; i++) 
    {
        node_ids[i] = i;
    }
    //create a random graph where each node will be connected to ateleast 3 nodes.
    vector<vector<pair<int,double>>> graph(n);
    srand(static_cast<unsigned int>(std::time(nullptr)));
    //start from the node 0
    for (int i = 0; i < n; i++) 
    {
        int j=0;
        random_shuffle(node_ids.begin(),node_ids.end());
        vector<int>peers;
        for(auto it: graph[i])
        {
        	peers.push_back(it.first);
        }
        //if the node has less than 3 peers randomly get a peer and establish a connection with the random pij as the latency to that connection.
        while (graph[i].size() < 3 && j<n) 
        {
            if (node_ids[j] != i && find(peers.begin(), peers.end(), node_ids[j])==peers.end()) 
            {
            	double pij=(double)(10+rand()%491)/1000;
            	graph[i].push_back({node_ids[j],pij});
            	graph[node_ids[j]].push_back({i,pij});
            }
            j++;
        }
    }
    //checking whether the created graph is connected or not using dfs.
    vector<int> visited(n,0);
    int flag=1;
    visited[0]=1;
    //call dfs from node 0, and store the visited array.
    dfs(0,graph,visited);
    //if any of the node is not visited after dfs, it means the graph is not connected.
    for(int i=0;i<n;i++)
    {
    	if(visited[i]==0)
    	{
    		flag=0;
    		break;
    	}
    }
    //if there is no unvisited node push the peers of each node to the network.
    if(flag)
    {
    	for(int i=0;i<n;i++)
    	{
    		network[i].Peers=graph[i];
    	}
    }
    //if there is any unvisited node create the graph from the start.
    else
    {
    	graph_creation(n);
    }
}
//function to find a block in blockchain starting from the given block using dfs.
Block* findBlockByIdDFS(Block *root, int targetId) {
    if (root == NULL) {
        return NULL;
    }

    if (root->id == targetId) {
        return root;
    }

    for (Block& pointedBy : root->pointed_by) 
    {
    	Block *result=findBlockByIdDFS(&pointedBy, targetId);
    	if(result)
    		return result;
    }

    return NULL;  // Block with targetId not found in the subtree 
}
//function to print the blockchain using dfs. 
void blockchain_traversal(Block * root)
{
	if(root==NULL) return;
	cout<<"Block id is"<<root->id<<" "<<root->time_added<<" parent is "<<root->parent_id <<endl;
	for(auto it:root->pointed_by)
	{
		blockchain_traversal(&it);
	}
}
//function to print the blockchain using bfs traversal.
void blockchain_traversal_bfs(Block* root,int node) {
    if (root == nullptr) return;

    std::queue<Block*> bfsQueue;
    bfsQueue.push(root);
    int t=1;
    Block* currentBlock;
    
    while (!bfsQueue.empty()) {
        currentBlock = bfsQueue.front();
        bfsQueue.pop();
        //cout<<endl;
        cout<<node<<","<<currentBlock->id<<","<<currentBlock->parent_id<<","<<currentBlock->time_added<<","<<currentBlock->length_of_chain<<","<<currentBlock->created_by<<endl;
        //std::cout <<t<< "Block id is " << currentBlock->id << " " << currentBlock->time_added << "created by "<<currentBlock->created_by<<" parent is " << currentBlock->parent_id <<" length_of_chain is"<<currentBlock->length_of_chain<< std::endl;
        t++;
        //cout<<endl;
        for (Block& neighbor : currentBlock->pointed_by) {
            bfsQueue.push(&neighbor);
        }
    }

}

int main(int argc,char *argv[])
{
	srand(static_cast<unsigned int>(std::time(nullptr))); 
	int z0=atoi(argv[1]);  //percentage of slow nodes
	int z1=atoi(argv[2]);  //percentage of low cpu nodes
	cout<<"Enter no.of nodes in the network";
	cin>>n;
	int last_block_id;
	transactions_gen_by_each_node.resize(n,0);
	blocks_gen_by_each_node.resize(n,0);
	blocks_added_by_each_node.resize(n,0);
	//creating a n size array of Node structures.
	network = new Node[n];
	//intialize the network and create the graph.
	initialization(n,z0,z1);
	graph_creation(n);

	//creating the genesis block which will be in each node at the start.
	//id of genesis block is 0 , parent id of genesis block is -1.
	//balances of each node in genesis block are 100.
	//length of the chain till genesis block is 1.
	Block genesis_block;
    genesis_block.id = 0;
    genesis_block.parent_id = -1;  // Assuming -1 as a special value for the genesis block
    genesis_block.time_added = 0;
    genesis_block.length_of_chain = 1;
    vector<int> balances(n,100);
    genesis_block.balance_nodes_at_this_block=balances;
    //adding the genesis block to each node in the network.
	for(int i=0;i<n;i++)
	{
		network[i].genesis_block=genesis_block;
		

	}
	//get exponential time with TX_MEAN and push a trasaction_generation action into event queue with that time.
	int time=get_exponential_dist_numb(1.0);
	pq.push(Event("txn_gen",0,0,time));
	//Every node will start mining the block at some random time with mean 100 sec. So push actions saying that they need to start first_block_gen at the choosen time.
	for(int i=0;i<n;i++)
	{
		int time=get_exponential_dist_numb(100);
		pq.push(Event("first_block_gen",0,i,time));
	}
	//loop until coin_base_transaction_id<100 i.e until 100 blocks are added to the block chain.
	while(!pq.empty())
	{
		auto action=pq.top();  //get the action with least time from the event queue and pop it.
		pq.pop();

		// if the action is firt block generation
		if(action.operation.compare("first_block_gen")==0)
		{
			//create a block with parent block as the block on which this node currently working on i.e intially the node will be pointing to genesis block as this is the first block node is creating
			Block new_block;
			new_block.balance_nodes_at_this_block=genesis_block.balance_nodes_at_this_block;

			new_block.created_by=action.src;
			new_block.id=global_block_id++;
			new_block.parent_id=network[action.src].working_on;
			int added=0;
			//pick the transactions which the node has seen till the block generation time and push them into the block.
			for(auto it:network[action.src].transactions)
			{
				if(added>100) break;  //if the no.of transactions pushed into block exceeds 100 stop adding further.
				Transaction txn=global_transactions[it.first];  
				if(it.second<100)
				{
					if(new_block.balance_nodes_at_this_block[txn.source]> txn.transfer_amount)  //if the transaction has transfer amount less than the node balance at this time i.e this transaction can be approved
					{

						new_block.balance_nodes_at_this_block[txn.source]-=txn.transfer_amount;   //change the balances of source and destination the transaction
						new_block.balance_nodes_at_this_block[txn.destination]+=txn.transfer_amount;
						new_block.transactions.push_back(it.first);   //push the transaction into block.
						new_block.transactions_till_now.push_back(it.first);  //push the transaction into transactions that are present in the chain till now.
						added++;
					}
				}
			}
			global_blocks[new_block.id]=new_block;   //push the block created into the global blocks map.
			int time=get_exponential_dist_numb((double)10/(network[action.src].hashing_power));  //get the interarrival time of the block from the exponential distribution time.
			pq.push(Event("block_gen",new_block.id,action.src,100+time));   //push the aciton block generation by this node with the time when it will crack POW.
			blocks_gen_by_each_node[action.src]++;
		}
		//if the action is transaction generation
		if(action.operation.compare("txn_gen")==0)
		{
			generate_transaction(action.time);   //generate random transaction at this time.
			int time=get_exponential_dist_numb(1.0);  //calculate the inter arrival time of transactions from the exponential distribution.
			pq.push(Event("txn_gen",0,0,action.time+time));  //push a transaction generation action with the new time.

		}
		//if the aciton is transaction receive , add the transaction to the node and broadcast it to its peers.
		else if(action.operation.compare("txn_recv")==0)
		{
			broadcast_transaction(action.id, action.src, action.time);

		}
		//if the action is block receive 
		else if(action.operation.compare("block_recv")==0)
		{
			//get the block with the block id in event.
			Block blk=global_blocks[action.id];
			//if the node has already seen this block continue 
			if(network[action.src].block_ids_seen.find(blk.id)!=network[action.src].block_ids_seen.end()) continue;
			//if the node has not seen this block before.
			//get the pointer to the parent to which this block is pointing it to
			Block *parent=findBlockByIdDFS(&network[action.src].genesis_block,blk.parent_id);
			if(parent)
			{
				//validating the block
				vector<int> balances=parent->balance_nodes_at_this_block;
				for(auto it:blk.transactions)
				{
					Transaction txn=global_transactions[it];
					balances[txn.source]-=txn.transfer_amount;
					balances[txn.destination]+=txn.transfer_amount;
				}
				int i=0;
				for( i=0;i<n;i++)
				{
					if(balances[i]!=blk.balance_nodes_at_this_block[i]) break;
				}
				//if the block is not valid i.e transactions in the block are not valid ignore this block.
				if(i!=n) 
				{
					continue;
				}
				//push this block into the block chain of this node.
				blk.time_added=action.time;
				blk.length_of_chain=parent->length_of_chain+1;
				parent->pointed_by.push_back(blk);
				//broadcast this block to it's peers.
				broadcast_block(action.id,action.src,action.time);
				//find the block this node is working on.
				Block * working_on_block=findBlockByIdDFS(&network[action.src].genesis_block, network[action.src].working_on);
				if(working_on_block)
				{
					//if the new block has created a longest chain than the chain which the node is currently working on.
					if(working_on_block->length_of_chain<blk.length_of_chain)
					{
						//change the mining  point of this node to the newly created longest chain.
						network[action.src].working_on=blk.id;

						//create a new block
						int added=0;
						if(coin_base_transaction_id==100) continue;
						Block new_block;
						new_block.transactions_till_now=blk.transactions_till_now;
						new_block.balance_nodes_at_this_block=blk.balance_nodes_at_this_block;
						for(auto it:network[action.src].transactions)
						{
							if(added>100) break;
							if(network[action.src].transactions[it.first]<action.time)
							{	
								Transaction txn=global_transactions[it.first];
								//if the transaction is not there in the chain in which we are adding this block.
								if(find(new_block.transactions_till_now.begin(),new_block.transactions_till_now.end(),it.first)==new_block.transactions_till_now.end())
								{
									//if the transaction is valid, then add this transaction to the block, and update balances of the source and destinatino of the transaction.
									if(new_block.balance_nodes_at_this_block[txn.source]>txn.transfer_amount)
									{
										new_block.balance_nodes_at_this_block[txn.source]-=txn.transfer_amount;
										new_block.balance_nodes_at_this_block[txn.destination]+=txn.transfer_amount;
										new_block.transactions.push_back(it.first);
										new_block.transactions_till_now.push_back(it.first);
										added++;
									}
								}
								
							}
						}
						
						new_block.id=global_block_id++;
						new_block.parent_id=network[action.src].working_on;   //the parent of this new block will be the block on which the node is currently working on.
						new_block.created_by=action.src;
						global_blocks[new_block.id]=new_block;  //add the block to the global blocks.
						int inter_arrival_time=get_exponential_dist_numb((double)10/ (network[action.src].hashing_power));
						//get the inter block arrival time and add the action block generation by this node at the new time.
						blocks_gen_by_each_node[action.src]++;
						pq.push(Event("block_gen",new_block.id,action.src,action.time+inter_arrival_time));
					}
				}	
				
			}
			
		}
		//if the action we got is block generation.
		else if(action.operation.compare("block_gen")==0)
		{
			//if the parent of this block is same as the block on which the node is currently mining on i.e new longest chain has not created till now.
			if(network[action.src].working_on==global_blocks[action.id].parent_id)
			{
				//get the block with that block id.
				Block new_block=global_blocks[action.id];
				//find the pointer to the parent block to which this block should be attached to.
				Block *parent=findBlockByIdDFS(&network[action.src].genesis_block,new_block.parent_id);
				if(parent)
				{
					blocks_added_by_each_node[action.src]++;
					//create a coin base transaction with 50 coins, and insert it into the start of the transactions of this block.
					Transaction txn;
					txn.id=coin_base_transaction_id;
					txn.source=action.src;
					txn.transfer_amount=50;
					txn.destination=action.src;
					string message= "TXN "+to_string(coin_base_transaction_id-1)+": User: "+to_string(action.src)+" GOT "+to_string(50)+" by mining";
					txn.message=message;
					//adding the coin base transaction to the global transactions.
					global_transactions[coin_base_transaction_id]=txn;
					new_block.transactions.insert(new_block.transactions.begin(),coin_base_transaction_id);  //inserting the coin base at the start of transactions of the block.
					global_blocks[action.id]=new_block;
					new_block.time_added=action.time;
					new_block.length_of_chain=parent->length_of_chain+1;  //length of the chain till this block is length of parent +1.

					parent->pointed_by.push_back(new_block);  //add this block to the parent childs.
					if(coin_base_transaction_id==99)
					{
						last_block_id=action.id;
						priority_queue< Event, vector<Event> ,compare_time> empty_queue;
						swap(pq,empty_queue);
					} 
					coin_base_transaction_id++;
					broadcast_block(action.id,action.src,action.time); //broadcast the block to it's peers.
				}
				
			}
			//if the block is not pointing to the block on which the node is currently working on,i.e we have got some new longest lenth chain.
			else{}
		}
	}

	//opening a file to push the output to visualize the blockchain.
	std::ofstream outputFile("output.csv", std::ios::trunc); 
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output.csv file!" << std::endl;
        return 0;
    }

    // Redirect cout to the output file
    std::streambuf* coutBuffer = std::cout.rdbuf();
    std::cout.rdbuf(outputFile.rdbuf());

    //print the blockchains maintained by all nodes.
	for(int id=0;id<n;id++)
	{
		cout<<"Node,Block ID,parent ID,Time seen,Length of chain,Created by"<<endl;
		blockchain_traversal_bfs(&network[id].genesis_block,id);
	}

	// Restore cout to the original buffer
    std::cout.rdbuf(coutBuffer);

    // Close the file
    outputFile.close();
    std::ofstream outputFile1("ex.csv", std::ios::trunc); 
    if (!outputFile1.is_open()) {
        std::cerr << "Error opening output.csv file!" << std::endl;
        return 0;
    }

    // Redirect cout to the output file
    std::streambuf* coutBuffer1 = std::cout.rdbuf();
    std::cout.rdbuf(outputFile1.rdbuf());
    cout<<"Node,type_node,type_cpu,txn_gen,blcK_gen,blck_added"<<endl;
    for(int i=0;i<n;i++)
    {
    	cout<<i<<","<<network[i].type_of_node<<","<<network[i].type_of_cpu<<","<<transactions_gen_by_each_node[i]<<","<<blocks_gen_by_each_node[i]<<","<<blocks_added_by_each_node[i]<<endl;
    }
    cout<<endl;
    std::cout.rdbuf(coutBuffer1);

    // Close the file
    outputFile1.close();
	
	
}
