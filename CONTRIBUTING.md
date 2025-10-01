# Contributing to TrackTx-NF

Thank you for your interest in contributing to TrackTx-NF! This document provides guidelines for contributing to the pipeline.

## 🚀 Quick Start for Contributors

### 1. Fork and Clone
```bash
git clone https://github.com/YOUR_USERNAME/TrackTx.git
cd TrackTx
```

### 2. Set Up Development Environment
```bash
# Install Docker (recommended)
# https://docs.docker.com/get-docker/

# Or install Nextflow and dependencies
# See README.md for full installation instructions
```

### 3. Test Your Changes
```bash
# Run tests with example data
./run_pipeline.sh --dry-run

# Test with small dataset
./run_pipeline.sh --samplesheet test_samplesheet.csv
```

## 📋 Contribution Guidelines

### Types of Contributions
- 🐛 **Bug fixes** - Fix issues in existing functionality
- ✨ **New features** - Add new analysis capabilities
- 📖 **Documentation** - Improve guides, README, code comments
- 🧪 **Tests** - Add test cases and validation
- 🔧 **Performance** - Optimize speed or memory usage

### Development Workflow

1. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   # or
   git checkout -b bugfix/issue-description
   ```

2. **Make your changes**
   - Follow existing code style
   - Add tests if applicable
   - Update documentation

3. **Test thoroughly**
   ```bash
   # Test with dry run
   ./run_pipeline.sh --dry-run
   
   # Test with small dataset
   ./run_pipeline.sh --samplesheet test_samplesheet.csv
   ```

4. **Commit with clear messages**
   ```bash
   git add .
   git commit -m "Add feature: brief description of what you added"
   git push origin feature/your-feature-name
   ```

5. **Open a Pull Request**
   - Provide clear description of changes
   - Link to any related issues
   - Include test results

## 🧪 Testing Guidelines

### Before Submitting
- [ ] Pipeline runs successfully with `--dry-run`
- [ ] All existing tests pass
- [ ] New functionality has test coverage
- [ ] Documentation is updated
- [ ] Code follows existing style

### Test Data
- Use small test datasets for development
- Don't commit large files (>10MB) to the repository
- Use `.gitignore` to exclude test data

## 📝 Code Style

### Python Files
- Follow PEP 8 style guide
- Use type hints where appropriate
- Add docstrings for functions
- Keep functions focused and small

### Nextflow Files
- Use consistent indentation (2 spaces)
- Add comments for complex logic
- Follow DSL2 best practices
- Use meaningful variable names

### Shell Scripts
- Use `set -euo pipefail`
- Quote variables: `"$variable"`
- Use `[[ ]]` for conditionals
- Add error handling

## 🐛 Reporting Bugs

### Before Reporting
1. Check existing [Issues](https://github.com/YOUR_USERNAME/TrackTx/issues)
2. Update to latest version
3. Try with minimal example

### Bug Report Template
```markdown
**Describe the bug**
A clear description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Run command '...'
2. See error

**Expected behavior**
What you expected to happen.

**Log files**
Include relevant log files (remove sensitive data).

**Environment:**
- OS: [e.g. Ubuntu 20.04]
- Nextflow version: [e.g. 25.04.7]
- Docker version: [e.g. 24.0.0]

**Additional context**
Any other context about the problem.
```

## ✨ Feature Requests

### Before Requesting
1. Check existing [Discussions](https://github.com/YOUR_USERNAME/TrackTx/discussions)
2. Search closed issues for similar requests
3. Consider if it fits the pipeline's scope

### Feature Request Template
```markdown
**Is your feature request related to a problem?**
A clear description of what the problem is.

**Describe the solution you'd like**
A clear description of what you want to happen.

**Describe alternatives you've considered**
Other solutions you've thought about.

**Additional context**
Any other context about the feature request.
```

## 📚 Documentation

### Types of Documentation
- **README.md** - Main documentation
- **docs/** - Detailed guides
- **Code comments** - Inline documentation
- **Examples** - Usage examples

### Documentation Standards
- Use clear, concise language
- Include code examples
- Add screenshots for UI elements
- Keep up-to-date with code changes

## 🔧 Development Setup

### Recommended Tools
- **IDE**: VS Code with Nextflow extension
- **Version Control**: Git with conventional commits
- **Testing**: Docker for consistent environments
- **Documentation**: Markdown with good formatting

### Environment Variables
```bash
# For development
export NXF_DEBUG=2
export NXF_ANSI_LOG=false

# For testing
export NXF_TEST_HOME=/tmp/nf-test
```

## 🤝 Code of Conduct

### Our Pledge
We are committed to providing a welcoming and inclusive experience for everyone, regardless of background, identity, or experience level.

### Expected Behavior
- Be respectful and inclusive
- Focus on constructive feedback
- Help others learn and grow
- Follow community guidelines

### Reporting Issues
If you experience or witness unacceptable behavior, please report it to the maintainers.

## 📞 Getting Help

### Community Support
- 💬 [GitHub Discussions](https://github.com/YOUR_USERNAME/TrackTx/discussions)
- 🐛 [GitHub Issues](https://github.com/YOUR_USERNAME/TrackTx/issues)
- 📖 [Documentation](https://github.com/YOUR_USERNAME/TrackTx/wiki)

### Maintainer Contact
- Direct message maintainers for urgent issues
- Use GitHub issues for bug reports
- Use discussions for questions and ideas

## 🏆 Recognition

Contributors will be recognized in:
- CONTRIBUTORS.md file
- Release notes
- Project documentation
- GitHub contributor graph

## 📄 License

By contributing to TrackTx-NF, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to TrackTx-NF! 🧬✨
