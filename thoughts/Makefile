CXX=g++
CPPFILES=main.cc ZielonkaTree.cc
CXXFLAGS=-Wall -Wextra
OBJDIR=./obj
OBJFILES=$(patsubst %.cpp,$(OBJDIR)/%.o, $(CPPFILES))
BINARY=bin

all: $(BINARY)

$(BINARY): $(OBJFILES)
	@echo "Creating binary.."
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $^

run: $(BINARY)
	@./$(BINARY)

clean:
	@echo "Removing object files.."
	rm -rf $(OBJDIR) 
	@echo "Removing binary.."
	rm -f $(BINARY)
